# download and install all julia dependencies required for the analysis 
julia src_jl/instantiate.jl
# parse wall source logfiles and save to hdf5
julia --project=./src_jl src_jl/parse.jl --ensemble_metadata ./metadata/ensembles.csv --output_hdf5 ./data_assets/wall_correlators.hdf5
# parse and analyse wilson flow/topology measurements 
echo -n "Parse and analyse wilson flow logs..."
python3 src_py/package_flows_multirep.py --h5_filename data_assets/topology.hdf5 --ensemble M1 ./raw_data/topology/Lt48Ls20beta6.5mf0.71mas1.01FUN/out/out_flow
python3 src_py/package_flows_multirep.py --h5_filename data_assets/topology.hdf5 --ensemble M2 ./raw_data/topology/Lt64Ls20beta6.5mf0.71mas1.01FUN/out/out_flow
python3 src_py/package_flows_multirep.py --h5_filename data_assets/topology.hdf5 --ensemble M3 ./raw_data/topology/Lt96Ls20beta6.5mf0.71mas1.01FUN/out/out_flow
python3 src_py/package_flows_multirep.py --h5_filename data_assets/topology.hdf5 --ensemble M4 ./raw_data/topology/Lt64Ls20beta6.5mf0.70mas1.01FUN/out/out_flow
python3 src_py/package_flows_multirep.py --h5_filename data_assets/topology.hdf5 --ensemble M5 ./raw_data/topology/Lt64Ls32beta6.5mf0.72mas1.01FUN/out/out_flow
echo " ...done!"
# determine autocorrelation times
echo -n "Analyse wilson flow logs..."
julia --project=./src_jl src_jl/autocorrelation.jl --wilson_flow_hdf5 data_assets/topology.hdf5 --wall_correlators_hdf5 data_assets/wall_correlators.hdf5 --output_hdf5 data_assets/autocor.hdf5 --output_tex assets/ensembles.tex --plot_dir data_assets/autocorrelation_plots 
echo " ...done!"
# perform individual fits for every ensemble and channel
echo -n "Fit wall source correlators..."
mkdir -p data_assets/wall_fits
python3 src_py/mass_wall.py --ensemble_name M1/FUN --plateau_start 16 --plateau_end 24  --output_file_mean data_assets/wall_fits/M1FUN_ps.csv --channel ps data_assets/wall_correlators.hdf5
python3 src_py/mass_wall.py --ensemble_name M1/FUN --plateau_start 16 --plateau_end 24  --output_file_mean data_assets/wall_fits/M1FUN_v.csv --channel v  data_assets/wall_correlators.hdf5
python3 src_py/mass_wall.py --ensemble_name M1/AS  --plateau_start 16 --plateau_end 24  --output_file_mean data_assets/wall_fits/M1AS_ps.csv  --channel ps data_assets/wall_correlators.hdf5
python3 src_py/mass_wall.py --ensemble_name M1/AS  --plateau_start 16 --plateau_end 24  --output_file_mean data_assets/wall_fits/M1AS_v.csv  --channel v  data_assets/wall_correlators.hdf5
python3 src_py/mass_wall.py --ensemble_name M2/FUN --plateau_start 22 --plateau_end 30  --output_file_mean data_assets/wall_fits/M2FUN_ps.csv --channel ps data_assets/wall_correlators.hdf5
python3 src_py/mass_wall.py --ensemble_name M2/FUN --plateau_start 22 --plateau_end 30  --output_file_mean data_assets/wall_fits/M2FUN_v.csv  --channel v data_assets/wall_correlators.hdf5
python3 src_py/mass_wall.py --ensemble_name M2/AS  --plateau_start 22 --plateau_end 30  --output_file_mean data_assets/wall_fits/M2AS_ps.csv  --channel ps data_assets/wall_correlators.hdf5
python3 src_py/mass_wall.py --ensemble_name M2/AS  --plateau_start 22 --plateau_end 30  --output_file_mean data_assets/wall_fits/M2AS_v.csv   --channel v data_assets/wall_correlators.hdf5
python3 src_py/mass_wall.py --ensemble_name M3/FUN --plateau_start 20 --plateau_end 30  --output_file_mean data_assets/wall_fits/M3FUN_v.csv --channel v  data_assets/wall_correlators.hdf5
python3 src_py/mass_wall.py --ensemble_name M3/FUN --plateau_start 20 --plateau_end 30  --output_file_mean data_assets/wall_fits/M3FUN_ps.csv --channel ps  data_assets/wall_correlators.hdf5
python3 src_py/mass_wall.py --ensemble_name M3/AS  --plateau_start 20 --plateau_end 30  --output_file_mean data_assets/wall_fits/M3AS_v.csv  --channel v  data_assets/wall_correlators.hdf5
python3 src_py/mass_wall.py --ensemble_name M3/AS  --plateau_start 20 --plateau_end 30  --output_file_mean data_assets/wall_fits/M3AS_ps.csv  --channel ps  data_assets/wall_correlators.hdf5
python3 src_py/mass_wall.py --ensemble_name M4/FUN --plateau_start 22 --plateau_end 32  --output_file_mean data_assets/wall_fits/M4FUN_ps.csv --channel ps data_assets/wall_correlators.hdf5
python3 src_py/mass_wall.py --ensemble_name M4/FUN --plateau_start 22 --plateau_end 32  --output_file_mean data_assets/wall_fits/M4FUN_v.csv --channel v  data_assets/wall_correlators.hdf5
python3 src_py/mass_wall.py --ensemble_name M4/AS  --plateau_start 20 --plateau_end 27  --output_file_mean data_assets/wall_fits/M4AS_ps.csv  --channel ps data_assets/wall_correlators.hdf5
python3 src_py/mass_wall.py --ensemble_name M4/AS  --plateau_start 20 --plateau_end 27  --output_file_mean data_assets/wall_fits/M4AS_v.csv  --channel v  data_assets/wall_correlators.hdf5
python3 src_py/mass_wall.py --ensemble_name M5/FUN --plateau_start 22 --plateau_end 30  --output_file_mean data_assets/wall_fits/M5FUN_ps.csv --channel ps data_assets/wall_correlators.hdf5
python3 src_py/mass_wall.py --ensemble_name M5/FUN --plateau_start 22 --plateau_end 30  --output_file_mean data_assets/wall_fits/M5FUN_v.csv  --channel v data_assets/wall_correlators.hdf5
python3 src_py/mass_wall.py --ensemble_name M5/AS  --plateau_start 22 --plateau_end 30  --output_file_mean data_assets/wall_fits/M5AS_ps.csv  --channel ps data_assets/wall_correlators.hdf5
python3 src_py/mass_wall.py --ensemble_name M5/AS  --plateau_start 22 --plateau_end 30  --output_file_mean data_assets/wall_fits/M5AS_v.csv   --channel v data_assets/wall_correlators.hdf5
echo " ...done!"
# create plots and table
echo -n "Create comparison plot and tables for decay constants..."
julia --project=./src_jl src_jl/compare.jl --wall_correlators_h5 data_assets/wall_correlators.hdf5 --wall_fits data_assets/wall_fits --smeared_results external_data/smeared --decay_output_csv data_assets/comparison_table.csv --decay_output_tex assets/local_smeared_decay_constants.tex --decay_output_pdf assets/wall_comparison.pdf
echo " ...done!"
echo "All done."