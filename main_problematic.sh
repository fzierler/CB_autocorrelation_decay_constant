python3 src_py/mass_wall.py --ensemble_name M2FUN --plateau_start 15 --plateau_end 25  --output_file_mean data_assets/M2FUN_p.csv --channel ps data_assets/test.hdf5
python3 src_py/mass_wall.py --ensemble_name M2AS  --plateau_start 15 --plateau_end 25  --output_file_mean data_assets/M2AS_p.csv  --channel ps data_assets/test.hdf5
python3 src_py/mass_wall.py --ensemble_name M3FUN --plateau_start 20 --plateau_end 40  --output_file_mean data_assets/M3FUN_v.csv --channel v  data_assets/test.hdf5
python3 src_py/mass_wall.py --ensemble_name M3AS  --plateau_start 20 --plateau_end 40  --output_file_mean data_assets/M3AS_v.csv  --channel v  data_assets/test.hdf5
python3 src_py/mass_wall.py --ensemble_name M5FUN --plateau_start 20 --plateau_end 27  --output_file_mean data_assets/M5FUN_p.csv --channel ps data_assets/test.hdf5
python3 src_py/mass_wall.py --ensemble_name M5AS  --plateau_start 20 --plateau_end 27  --output_file_mean data_assets/M5AS_p.csv  --channel ps data_assets/test.hdf5