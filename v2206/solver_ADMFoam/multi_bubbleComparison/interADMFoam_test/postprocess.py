import re
import pandas as pd

# Input/output files
input_file = "log.interADMFoam"
output_file = "interADMFoam_test_v4.csv"

# Regex patterns
time_pattern = re.compile(r"^\s*Time\s*=\s*([0-9.eE+-]+)\s*$")
gaseous_pattern = re.compile(
    r">>>\s*(\w+)\s+concentration\s*=\s*([0-9.eE+-]+),\s*generation rate\s*=\s*([0-9.eE+-]+)"
)
soluable_pattern = re.compile(
    # r">>>\s*(\w+)\s*=\s*([0-9.eE+-]+)$"
    r">>>\s*(\w+)\s*=\s*([0-9.eE+-]+)\s*$"
)

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

# Reorder columns
expected_cols = ["Time", "Sh2Ave", "Sch4Ave", "Sco2Ave", "Gh2", "RGh2", "Gch4", "RGch4", "Gco2", "RGco2"]
df = df.reindex(columns=expected_cols)

# Save to CSV
# df.to_csv(output_file, index=False)
df_sampled = df.iloc[::100]
# Ensure last row is included
if df.index[-1] != df_sampled.index[-1]:
    df_sampled = pd.concat([df_sampled, df.iloc[[-1]]])

df_sampled.to_csv(output_file, index=False)

print(f"Extracted {len(df)} time steps and saved to {output_file}")
