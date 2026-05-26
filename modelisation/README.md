Input data :
    * xp_input_NOT_TOUCHED : Row files shared by Lucien
    * xp_input_MODIFIED : Uniformised files, used in generate_data.jl to generate xp_input_20°
    * xp_input_20° : flow cytometer data for T=20°C, one file per (replicate x cycle)

Codes :
    * models.jl : Model definition
        * S
        * I
        * R
    * 030426_code_optimise_simulate_[fit-plot-model].jl : Run one model with fixed or fitted parameters
    * 040426_code_test_parameters_variability.jl : Test the variability of parameters


Communication :
    * README.md : presents the current architecture
    * README_day_to_day.md : what is done each day
    * photos/ : photos of board after meetings