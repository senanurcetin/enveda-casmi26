"""Fresh supervised fingerprint fitting shared by the external study runner."""
from pathlib import Path

def fit_fingerprints(store, labels, splits, protocol, protocol_path, artifact_dir, device):
    import numpy as np
    import torch
    from rdkit import Chem, DataStructs
    from rdkit.Chem import rdFingerprintGenerator
    from casmi.channels.fpmodel import FPModel, FPModelConfig, save
    from casmi.channels.fpinfer import FP_CLEAN, MAX_PEAKS, ADDUCT_ID
    from casmi.spectra.clean import clean_peaks_fast
    from casmi.validate.scaffold_study import digest
    from casmi.validate.source_manifest import sha256_file
    meta = store.meta
    key_of_row = meta.train_ik14.astype(str).to_numpy()
    seed = protocol['seed']
    artifact_dir = Path(artifact_dir)
    # Both D3 models train on selected training identities only, on opposite folds.
    train_lookup = {row["key"]: row for row in splits["train"]}
    fp_eligible = (meta.test_adduct.to_numpy() & (meta.n_peaks.to_numpy() >= 5)
                   & (np.abs(meta.ppm.to_numpy()) <= 20))
    training_rows = {}
    for index in np.flatnonzero(fp_eligible & np.isin(key_of_row, list(train_lookup))):
        training_rows.setdefault(key_of_row[index], []).append(int(index))
    target_gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    fp_models, fp_lineage = {}, []
    settings = protocol["fp_model"]
    (artifact_dir / "models").mkdir(exist_ok=True)
    for fold in (0, 1):
        rng = np.random.default_rng(seed + fold)
        selected_rows, selected_keys = [], []
        for key, item in sorted(train_lookup.items()):
            if item["fold"] != fold:
                continue
            available = np.asarray(training_rows.get(key, []), dtype=np.int64)
            cap = protocol["max_training_spectra_per_structure"]
            if len(available) > cap:
                available = np.sort(rng.choice(available, cap, replace=False))
            selected_rows.extend(available.tolist())
            selected_keys.extend([key] * len(available))
        if len(selected_rows) < 100:
            raise ValueError("Insufficient fingerprint training spectra in declared fold")
        if not set(selected_keys) <= {r["key"] for r in splits["train"]}:
            raise ValueError("Fingerprint training crossed the declared training partition")
        n = len(selected_rows)
        mz = np.zeros((n, MAX_PEAKS), np.float32)
        it = np.zeros_like(mz)
        npk = np.zeros(n, np.int64)
        unique_keys = sorted(set(selected_keys))
        targets = {}
        for key in unique_keys:
            mol = Chem.MolFromSmiles(labels[key])
            bits = np.zeros(2048, np.uint8)
            DataStructs.ConvertToNumpyArray(target_gen.GetFingerprint(mol), bits)
            if not bits.any():
                raise ValueError("Invalid fingerprint training target")
            targets[key] = bits
        for j, index in enumerate(selected_rows):
            cm, ci = clean_peaks_fast(*store.peaks(index), float(meta.precursor_mz.iat[index]), FP_CLEAN)
            mz[j, :len(cm)], it[j, :len(ci)], npk[j] = cm, ci, len(cm)
        side = meta.iloc[selected_rows]
        inputs = {"mz": mz, "it": it, "npk": npk, "prec": side.precursor_mz.to_numpy(np.float32),
                  "adduct": side.adduct.map(ADDUCT_ID).to_numpy(np.int64),
                  "ce": side.ce_ev.to_numpy(np.float32),
                  "positive": (side.ionization_mode == "positive").to_numpy(np.float32)}
        y = np.asarray([targets[key] for key in selected_keys], dtype=np.float32)
        torch.manual_seed(seed + fold)
        model = FPModel(FPModelConfig(d_model=settings["d_model"], n_layers=settings["n_layers"],
                                      d_ff=settings["d_ff"])).to(device)
        optimizer = torch.optim.AdamW(model.parameters(), lr=settings["lr"], weight_decay=.01)
        losses = []
        model.train()
        for epoch in range(settings["epochs"]):
            total = 0.0
            order = rng.permutation(n)
            for start in range(0, n, settings["batch_size"]):
                ids = order[start:start + settings["batch_size"]]
                x = {name: torch.as_tensor(np.ascontiguousarray(value[ids]), device=device,
                     dtype=torch.long if name in ("npk", "adduct") else torch.float32) for name, value in inputs.items()}
                target = torch.as_tensor(y[ids], device=device)
                loss = torch.nn.functional.binary_cross_entropy_with_logits(model(**x), target)
                if not torch.isfinite(loss):
                    raise ValueError("Nonfinite fingerprint training loss")
                optimizer.zero_grad(set_to_none=True)
                loss.backward()
                torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
                optimizer.step()
                total += float(loss.detach().cpu()) * len(ids)
            losses.append(total / n)
            print(f"Trained fresh fingerprint fold {fold}, epoch {epoch + 1}/{settings['epochs']}", flush=True)
        model.eval()
        model_path = artifact_dir / "models" / f"external_fp_fold{fold}.pt"
        lineage = {"fold": fold, "n_spectra": n, "n_structures": len(unique_keys),
                   "training_rows_sha256": digest(selected_rows), "training_identities_sha256": digest(unique_keys),
                   "training_families_sha256": digest(sorted({train_lookup[k]["group"] for k in unique_keys})),
                   "dev_test_identities_used": 0, "epochs_fixed_before_execution": settings["epochs"],
                   "training_bce_by_epoch": losses, "initialization": "fresh seeded random; no archived weights"}
        save(model, model_path, {"lineage": lineage, "protocol_sha256": sha256_file(protocol_path)})
        fp_lineage.append({**lineage, "model_sha256": sha256_file(model_path)})
        fp_models[fold] = model
        del mz, it, npk, inputs, y

    return fp_models, fp_lineage
