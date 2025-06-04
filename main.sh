# download and install all julia dependencies required for the analysis 
julia src_jl/instantiate.jl
# parse wall source logfiles and save to hdf5
julia src_jl/parse.jl
# parse and analyse wilson flow/topology measurements 
python3 src_py/package_flows_multirep.py --h5_filename data_assets/topology.hdf5 --ensemble M1 ./raw_data/topology/Lt48Ls20beta6.5mf0.71mas1.01FUN/out/out_flow
python3 src_py/package_flows_multirep.py --h5_filename data_assets/topology.hdf5 --ensemble M2 ./raw_data/topology/Lt64Ls20beta6.5mf0.71mas1.01FUN/out/out_flow
python3 src_py/package_flows_multirep.py --h5_filename data_assets/topology.hdf5 --ensemble M3 ./raw_data/topology/Lt96Ls20beta6.5mf0.71mas1.01FUN/out/out_flow
python3 src_py/package_flows_multirep.py --h5_filename data_assets/topology.hdf5 --ensemble M4 ./raw_data/topology/Lt64Ls20beta6.5mf0.70mas1.01FUN/out/out_flow
python3 src_py/package_flows_multirep.py --h5_filename data_assets/topology.hdf5 --ensemble M5 ./raw_data/topology/Lt64Ls32beta6.5mf0.72mas1.01FUN/out/out_flow
# determine autocorrelation times
julia src_jl/autocorrelation.jl    
# perform individual fits for every ensemble and channel
python3 src_py/mass_wall.py --ensemble_name M1/FUN --plateau_start 16 --plateau_end 24  --output_file_mean data_assets/M1FUN_ps.csv --channel ps data_assets/wall_correlators.hdf5
python3 src_py/mass_wall.py --ensemble_name M1/FUN --plateau_start 16 --plateau_end 24  --output_file_mean data_assets/M1FUN_v.csv --channel v  data_assets/wall_correlators.hdf5
python3 src_py/mass_wall.py --ensemble_name M1/AS  --plateau_start 16 --plateau_end 24  --output_file_mean data_assets/M1AS_ps.csv  --channel ps data_assets/wall_correlators.hdf5
python3 src_py/mass_wall.py --ensemble_name M1/AS  --plateau_start 16 --plateau_end 24  --output_file_mean data_assets/M1AS_v.csv  --channel v  data_assets/wall_correlators.hdf5
python3 src_py/mass_wall.py --ensemble_name M2/FUN --plateau_start 22 --plateau_end 30  --output_file_mean data_assets/M2FUN_ps.csv --channel ps data_assets/wall_correlators.hdf5
python3 src_py/mass_wall.py --ensemble_name M2/FUN --plateau_start 22 --plateau_end 30  --output_file_mean data_assets/M2FUN_v.csv  --channel v data_assets/wall_correlators.hdf5
python3 src_py/mass_wall.py --ensemble_name M2/AS  --plateau_start 22 --plateau_end 30  --output_file_mean data_assets/M2AS_ps.csv  --channel ps data_assets/wall_correlators.hdf5
python3 src_py/mass_wall.py --ensemble_name M2/AS  --plateau_start 22 --plateau_end 30  --output_file_mean data_assets/M2AS_v.csv   --channel v data_assets/wall_correlators.hdf5
python3 src_py/mass_wall.py --ensemble_name M3/FUN --plateau_start 20 --plateau_end 30  --output_file_mean data_assets/M3FUN_v.csv --channel v  data_assets/wall_correlators.hdf5
python3 src_py/mass_wall.py --ensemble_name M3/FUN --plateau_start 20 --plateau_end 30  --output_file_mean data_assets/M3FUN_ps.csv --channel ps  data_assets/wall_correlators.hdf5
python3 src_py/mass_wall.py --ensemble_name M3/AS  --plateau_start 20 --plateau_end 30  --output_file_mean data_assets/M3AS_v.csv  --channel v  data_assets/wall_correlators.hdf5
python3 src_py/mass_wall.py --ensemble_name M3/AS  --plateau_start 20 --plateau_end 30  --output_file_mean data_assets/M3AS_ps.csv  --channel ps  data_assets/wall_correlators.hdf5
python3 src_py/mass_wall.py --ensemble_name M4/FUN --plateau_start 22 --plateau_end 32  --output_file_mean data_assets/M4FUN_ps.csv --channel ps data_assets/wall_correlators.hdf5
python3 src_py/mass_wall.py --ensemble_name M4/FUN --plateau_start 22 --plateau_end 32  --output_file_mean data_assets/M4FUN_v.csv --channel v  data_assets/wall_correlators.hdf5
python3 src_py/mass_wall.py --ensemble_name M4/AS  --plateau_start 20 --plateau_end 27  --output_file_mean data_assets/M4AS_ps.csv  --channel ps data_assets/wall_correlators.hdf5
python3 src_py/mass_wall.py --ensemble_name M4/AS  --plateau_start 20 --plateau_end 27  --output_file_mean data_assets/M4AS_v.csv  --channel v  data_assets/wall_correlators.hdf5
python3 src_py/mass_wall.py --ensemble_name M5/FUN --plateau_start 22 --plateau_end 30  --output_file_mean data_assets/M5FUN_ps.csv --channel ps data_assets/wall_correlators.hdf5
python3 src_py/mass_wall.py --ensemble_name M5/FUN --plateau_start 22 --plateau_end 30  --output_file_mean data_assets/M5FUN_v.csv  --channel v data_assets/wall_correlators.hdf5
python3 src_py/mass_wall.py --ensemble_name M5/AS  --plateau_start 22 --plateau_end 30  --output_file_mean data_assets/M5AS_ps.csv  --channel ps data_assets/wall_correlators.hdf5
python3 src_py/mass_wall.py --ensemble_name M5/AS  --plateau_start 22 --plateau_end 30  --output_file_mean data_assets/M5AS_v.csv   --channel v data_assets/wall_correlators.hdf5
## create plot and table
julia src_jl/compare.jl