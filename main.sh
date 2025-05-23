# parse logfiles and save to hdf5
julia parse.jl
# create plots of the effective mass for further reference
julia effective_mass.jl
# perform individual fits for every ensemble and channel
python3 src_py/mass_wall.py --ensemble_name M1FUN --plateau_start 16 --plateau_end 24  --output_file_mean data_assets/M1FUN_ps.csv --channel ps data_assets/test.hdf5
python3 src_py/mass_wall.py --ensemble_name M1FUN --plateau_start 16 --plateau_end 24  --output_file_mean data_assets/M1FUN_v.csv --channel v  data_assets/test.hdf5
python3 src_py/mass_wall.py --ensemble_name M1AS  --plateau_start 16 --plateau_end 24  --output_file_mean data_assets/M1AS_ps.csv  --channel ps data_assets/test.hdf5
python3 src_py/mass_wall.py --ensemble_name M1AS  --plateau_start 16 --plateau_end 24  --output_file_mean data_assets/M1AS_v.csv  --channel v  data_assets/test.hdf5
python3 src_py/mass_wall.py --ensemble_name M2FUN --plateau_start 22 --plateau_end 30  --output_file_mean data_assets/M2FUN_ps.csv --channel ps data_assets/test.hdf5
python3 src_py/mass_wall.py --ensemble_name M2FUN --plateau_start 22 --plateau_end 30  --output_file_mean data_assets/M2FUN_v.csv  --channel v data_assets/test.hdf5
python3 src_py/mass_wall.py --ensemble_name M2AS  --plateau_start 22 --plateau_end 30  --output_file_mean data_assets/M2AS_ps.csv  --channel ps data_assets/test.hdf5
python3 src_py/mass_wall.py --ensemble_name M2AS  --plateau_start 22 --plateau_end 30  --output_file_mean data_assets/M2AS_v.csv   --channel v data_assets/test.hdf5
python3 src_py/mass_wall.py --ensemble_name M3FUN --plateau_start 20 --plateau_end 30  --output_file_mean data_assets/M3FUN_v.csv --channel v  data_assets/test.hdf5
python3 src_py/mass_wall.py --ensemble_name M3FUN --plateau_start 20 --plateau_end 30  --output_file_mean data_assets/M3FUN_ps.csv --channel ps  data_assets/test.hdf5
python3 src_py/mass_wall.py --ensemble_name M3AS  --plateau_start 20 --plateau_end 30  --output_file_mean data_assets/M3AS_v.csv  --channel v  data_assets/test.hdf5
python3 src_py/mass_wall.py --ensemble_name M3AS  --plateau_start 20 --plateau_end 30  --output_file_mean data_assets/M3AS_ps.csv  --channel ps  data_assets/test.hdf5
python3 src_py/mass_wall.py --ensemble_name M4FUN --plateau_start 22 --plateau_end 32  --output_file_mean data_assets/M4FUN_ps.csv --channel ps data_assets/test.hdf5
python3 src_py/mass_wall.py --ensemble_name M4FUN --plateau_start 22 --plateau_end 32  --output_file_mean data_assets/M4FUN_v.csv --channel v  data_assets/test.hdf5
python3 src_py/mass_wall.py --ensemble_name M4AS  --plateau_start 20 --plateau_end 27  --output_file_mean data_assets/M4AS_ps.csv  --channel ps data_assets/test.hdf5
python3 src_py/mass_wall.py --ensemble_name M4AS  --plateau_start 20 --plateau_end 27  --output_file_mean data_assets/M4AS_v.csv  --channel v  data_assets/test.hdf5
python3 src_py/mass_wall.py --ensemble_name M5FUN --plateau_start 22 --plateau_end 30  --output_file_mean data_assets/M5FUN_ps.csv --channel ps data_assets/test.hdf5
python3 src_py/mass_wall.py --ensemble_name M5FUN --plateau_start 22 --plateau_end 30  --output_file_mean data_assets/M5FUN_v.csv  --channel v data_assets/test.hdf5
python3 src_py/mass_wall.py --ensemble_name M5AS  --plateau_start 22 --plateau_end 30  --output_file_mean data_assets/M5AS_ps.csv  --channel ps data_assets/test.hdf5
python3 src_py/mass_wall.py --ensemble_name M5AS  --plateau_start 22 --plateau_end 30  --output_file_mean data_assets/M5AS_v.csv   --channel v data_assets/test.hdf5
# create plot and table
julia compare.jl