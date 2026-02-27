import re
import pandas as pd

# Input/output files
input_file = "log.interADMFoam"
output_file = "batch_mixedNoBubbles.csv"

# Regex patterns
time_pattern = re.compile(r"^\s*Time\s*=\s*([0-9.eE+-]+)\s*$")
gaseous_pattern = re.compile(
    r">>>\s*(\w+)\s+concentration\s*=\s*([0-9.eE+-]+),\s*generation rate\s*=\s*([0-9.eE+-]+)"
)
soluable_pattern = re.compile(
    # r">>>\s*(\w+)\s*=\s*([0-9.eE+-]+)$"
    r">>>\s*(\w+)\s*=\s*([0-9.eE+-]+)\s*$"
)

# --- collect setup lines from the log (will be written as comments at top of CSV)
setup_keys = [
    "is benchmark case",
    "Qin_",
    "Vgas_",
    "Vliq_",
    "Vfrac_",
    "qGas_",
]
setup_lines = []
try:
    with open(input_file, "r") as _lf:
        for _l in _lf:
            if _l.startswith(">>>"):
                s = _l[3:].strip()
                for _k in setup_keys:
                    if s.startswith(_k):
                        # keep the original trimmed line (e.g. 'is benchmark case: 0')
                        setup_lines.append(s)
                        break
except Exception:
    # if the log isn't available or readable, leave setup_lines empty
    setup_lines = []

# --- override or add benchmark value from constant/admno1Properties
try:
    adm_path = "constant/admno1Properties"
    with open(adm_path, "r") as _adm:
        adm_text = _adm.read()
    bm_match = re.search(r"^\s*benchmark\s+([^;\s]+)", adm_text, re.MULTILINE)
    if bm_match:
        bm_val = bm_match.group(1).strip()
        # replace any existing 'is benchmark case' entry or add it at the top
        replaced = False
        for i, _s in enumerate(setup_lines):
            if _s.startswith("is benchmark case"):
                setup_lines[i] = f"Benchmark: {bm_val}"
                replaced = True
                break
        if not replaced:
            setup_lines.insert(0, f"Benchmark: {bm_val}")
except Exception:
    # ignore if admno1Properties missing or unreadable
    pass

# --- read alpha values from constant/transportProperties
try:
    tp_path = "constant/transportProperties"
    with open(tp_path, "r") as _tp:
        for _l in _tp:
            m = re.search(r"^\s*(alphaW|alphaI|alphaInterface)\s+([0-9.eE+-]+)", _l)
            if m:
                key = m.group(1)
                val = m.group(2)
                # only add once and keep order after benchmark
                entry = f"{key}: {val}"
                if entry not in setup_lines:
                    setup_lines.append(entry)
except Exception:
    # ignore missing transportProperties
    pass
    

# Storage
records = []
current_time = None
current_record = {}

with open(input_file, "r") as f:
    for line in f:
        # Detect time
        tmatch = time_pattern.search(line)
        if tmatch:
            # If we already collected a block, store it
            if current_record:
                records.append(current_record)
                current_record = {}
            current_time = float(tmatch.group(1))
            current_record["Time"] = current_time
            continue

        # Detect soluable
        smatch = soluable_pattern.search(line)
        if smatch and current_time is not None:
            species = smatch.group(1)   # e.g. Gh2
            conc = float(smatch.group(2))
            current_record[species] = conc
            continue

        # Detect gaseous
        dmatch = gaseous_pattern.search(line)
        if dmatch and current_time is not None:
            species = dmatch.group(1)   # e.g. Gh2
            conc = float(dmatch.group(2))
            rate = float(dmatch.group(3))
            current_record[species] = conc
            current_record[f"R{species}"] = rate

        

# Don’t forget last block
if current_record:
    records.append(current_record)

# Convert to DataFrame
df = pd.DataFrame(records)

# Reorder columns (read with 'Ave' suffix, write without)
read_cols = ["Time", "Sh2Ave", "Sch4Ave", "Sco2Ave", "Gh2", "RGh2", "Gch4", "RGch4", "Gco2", "RGco2"]
write_cols = ["Time", "Sh2", "Sch4", "Sco2", "Gh2", "RGh2", "Gch4", "RGch4", "Gco2", "RGco2"]
df = df.reindex(columns=read_cols)

# Rename columns for output
rename_map = {r: w for r, w in zip(read_cols, write_cols) if r != w}
df = df.rename(columns=rename_map)


# Save to CSV
# df.to_csv(output_file, index=False)
df_sampled = df.iloc[::100]
# Ensure last row is included
if df.index[-1] != df_sampled.index[-1]:
    df_sampled = pd.concat([df_sampled, df.iloc[[-1]]])

df_sampled.to_csv(output_file, index=False)

print(f"Extracted {len(df)} time steps and saved to {output_file}")

# --- prepend setup comments to the CSV output (if any setup lines were found)
if setup_lines:
    try:
        with open(output_file, "r") as _orig:
            _csv = _orig.read()
        with open(output_file, "w") as _out:
            _out.write("# Setup from log.interADMFoam\n")
            for _s in setup_lines:
                _out.write(f"# {_s}\n")
            _out.write("\n")
            _out.write(_csv)
    except Exception as _e:
        print(f"Warning: could not prepend setup header to {output_file}: {_e}")
