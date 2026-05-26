using XLSX
using DataFrames
using Dates
using CSV
using Base.Threads

# =========================
# Directories
# =========================
dir_input  = "modelisation/input/xp_input_MODIFIED"
dir_output = "modelisation/input/xp_input_all"
mkpath(dir_output)

# =========================
# Parameters
# =========================
datas      = ["host", "virus"]
conditions = ["host", "coevo", "virus"]
temps      = ["15", "20", "26"]
replicates = ["A", "B", "C"]
cycles     = [1, 2, 3, 4, 5, 6, 7]

# =========================
# Helpers
# =========================
compute_time_hours(dates) = (dates .- dates[1]) ./ Hour(1)

# =========================
# Main loop
# =========================
for data in datas, condition in conditions, temp in temps
    file_path = joinpath(dir_input, "data_$(data).xlsx")
    if !isfile(file_path)
        println("⚠️ File not found: $(file_path)")
        continue
    end

    xf = XLSX.readxlsx(file_path)
    sheet_name = "$(condition)_$(temp)"
    if !(sheet_name in XLSX.sheetnames(xf))
        println("⚠️ Sheet not found: $(sheet_name) in $(file_path)")
        continue
    end

    # Lecture unique du fichier
    df = DataFrame(XLSX.readtable(file_path, sheet_name))
    df.Date = DateTime.(df.Date)

    # Préparer toutes les combinaisons de cycle × replicate
    tasks = [(cycle, replicate) for cycle in cycles, replicate in replicates]

    # Parallélisation
    @threads for t in tasks
        cycle, replicate = t
        mask = (df.Cycle .== cycle) .& .!ismissing.(df.Cycle)

        out = DataFrame(
            time = compute_time_hours(df.Date)[mask],
            X    = df[mask, Symbol(replicate)]
        )

        output_file = joinpath(
            dir_output,
            "$(data)Data_$(condition)Condition_Temperature$(temp)_Replicate$(replicate)_Cycle$(cycle).csv"
        )
        CSV.write(output_file, out)
        println("✓ Written: $(basename(output_file))")
    end
end
