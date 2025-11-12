#!/usr/bin/env python
"""
Test program to compare fetch_station_product results with different datum values (MSL vs NAVD)
for NOAA station 8651370 (Wrightsville, NC)
"""

import sys
import os
import numpy as np
import pandas as pd
import datetime as dt

# Add the project root to the path so we can import modules
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from harvester.get_observations_stations import get_obs_stations
from utilities.utilities import utilities as utilities

def test_datum_comparison():
    """
    Test fetch_station_product with MSL and NAVD datums for station 8651370
    """
    # Initialize logging with the provided config file
    print("Initializing logging...")
    sys.stdout.flush()
    config_file = os.path.expanduser('~/floodwater_configs/suites/NCSCv2/gfs/NCSCv2.0_gfs/support_files/data_assimilation.yaml')
    print(f"Using config file: {config_file}")
    sys.stdout.flush()
    try:
        utilities.init_logging(subdir=None, config_file=config_file)
        print("✓ Logging initialized successfully")
        sys.stdout.flush()
    except Exception as e:
        print(f"✗ Error initializing logging: {e}")
        import traceback
        traceback.print_exc()
        sys.stdout.flush()
        return
    
    # Station ID for Wrightsville, NC
    station_id = '8658163'
    
    # Set up time range: 11/01/2025 00:00 UTC to 11/08/2025 00:00 UTC
    start_time = dt.datetime(2025, 11, 1, 0, 0, 0)
    end_time = dt.datetime(2025, 11, 8, 0, 0, 0)
    
    starttime_str = start_time.strftime('%Y-%m-%d %H:%M:%S')
    endtime_str = end_time.strftime('%Y-%m-%d %H:%M:%S')
    
    print(f"\n{'='*80}")
    print(f"Testing datum comparison for station {station_id} (Wrightsville, NC)")
    print(f"Time range: {starttime_str} to {endtime_str}")
    print(f"{'='*80}\n")
    sys.stdout.flush()
    
    # Create get_obs_stations instance with NOAAWEB source
    # Note: We need to pass a station_list_file, but we'll override it anyway
    print("Creating get_obs_stations instance...")
    sys.stdout.flush()
    try:
        # Use a minimal station file with just our test station
        # The file will be overridden anyway, but get_noaa_stations requires a valid file
        test_station_file = os.path.join(os.path.dirname(__file__), 'test_station_list.csv')
        if not os.path.exists(test_station_file):
            # Create it if it doesn't exist
            # Format: header row, then skip row (empty), then data rows
            # The code uses skiprows=[1], so we need: header, empty row, data
            with open(test_station_file, 'w') as f:
                f.write('serial_nr,stationid\n')
                f.write('\n')  # Empty row that will be skipped
                f.write('0,8658163\n')
        station_list_file = test_station_file
        
        obs = get_obs_stations(
            source='NOAAWEB',
            product='water_level',
            contrails_yamlname=None,
            knockout_dict=None,
            station_list_file=station_list_file  # Will be overridden
        )
        print(f"✓ get_obs_stations instance created successfully")
        sys.stdout.flush()
    except Exception as e:
        print(f"✗ Error creating get_obs_stations instance: {e}")
        import traceback
        traceback.print_exc()
        sys.stdout.flush()
        return
    
    # Override station list to only include our test station
    print(f"Overriding station list to include only {station_id}...")
    sys.stdout.flush()
    try:
        obs.override_station_IDs([station_id])
        print(f"✓ Station list set to: {obs.station_list}\n")
        sys.stdout.flush()
    except Exception as e:
        print(f"✗ Error overriding station list: {e}")
        import traceback
        traceback.print_exc()
        return
    
    # Test with MSL datum
    print("Fetching data with datum='MSL'...")
    try:
        data_msl, meta_msl = obs.fetch_station_product(
            (starttime_str, endtime_str),
            return_sample_min=0,
            interval=None,
            datum='MSL'
        )
        print(f"MSL fetch completed successfully")
        print(f"MSL data shape: {data_msl.shape if not isinstance(data_msl, float) else 'No data'}")
        if not isinstance(data_msl, float) and not data_msl.empty:
            print(f"MSL data columns: {data_msl.columns.tolist()}")
            print(f"MSL data index range: {data_msl.index.min()} to {data_msl.index.max()}")
            print(f"MSL sample values (first 5):")
            print(data_msl.head())
    except Exception as e:
        print(f"Error fetching MSL data: {e}")
        import traceback
        traceback.print_exc()
        data_msl = None
        meta_msl = None
    
    print("\n" + "-"*80 + "\n")
    
    # Test with NAVD datum
    print("Fetching data with datum='NAVD'...")
    try:
        data_navd, meta_navd = obs.fetch_station_product(
            (starttime_str, endtime_str),
            return_sample_min=0,
            interval=None,
            datum='NAVD'
        )
        print(f"NAVD fetch completed successfully")
        print(f"NAVD data shape: {data_navd.shape if not isinstance(data_navd, float) else 'No data'}")
        if not isinstance(data_navd, float) and not data_navd.empty:
            print(f"NAVD data columns: {data_navd.columns.tolist()}")
            print(f"NAVD data index range: {data_navd.index.min()} to {data_navd.index.max()}")
            print(f"NAVD sample values (first 5):")
            print(data_navd.head())
    except Exception as e:
        print(f"Error fetching NAVD data: {e}")
        import traceback
        traceback.print_exc()
        data_navd = None
        meta_navd = None
    
    print("\n" + "="*80 + "\n")
    
    # Compare the results
    print("COMPARISON RESULTS:")
    print("="*80)
    
    if data_msl is None or isinstance(data_msl, float):
        print("MSL data is not available for comparison")
    elif data_navd is None or isinstance(data_navd, float):
        print("NAVD data is not available for comparison")
    else:
        # Check if both have data
        if data_msl.empty:
            print("MSL data is empty")
        elif data_navd.empty:
            print("NAVD data is empty")
        else:
            # Align the dataframes on their time index
            # Find common time indices
            common_times = data_msl.index.intersection(data_navd.index)
            
            if len(common_times) == 0:
                print("No common time indices between MSL and NAVD data")
            else:
                print(f"\nNumber of common time points: {len(common_times)}")
                
                # Extract values for the station
                msl_values = data_msl.loc[common_times, station_id]
                navd_values = data_navd.loc[common_times, station_id]
                
                # Calculate the difference
                difference = navd_values - msl_values
                
                print(f"\nDifference statistics (NAVD - MSL):")
                print(f"  Mean difference: {difference.mean():.6f} meters")
                print(f"  Std deviation: {difference.std():.6f} meters")
                print(f"  Min difference: {difference.min():.6f} meters")
                print(f"  Max difference: {difference.max():.6f} meters")
                
                # Check if the difference is constant (datum offset)
                if difference.std() < 1e-6:
                    print(f"\n✓ The difference is constant: {difference.iloc[0]:.6f} meters")
                    print(f"  This indicates a fixed datum offset between NAVD and MSL")
                else:
                    print(f"\n⚠ The difference varies, which may indicate:")
                    print(f"  - Different data processing")
                    print(f"  - Different data sources")
                    print(f"  - Or time-varying datum corrections")
                
                # Show some sample comparisons
                print(f"\nSample comparison (first 10 common time points):")
                comparison_df = pd.DataFrame({
                    'Time': common_times[:10],
                    'MSL (m)': msl_values.iloc[:10].values,
                    'NAVD (m)': navd_values.iloc[:10].values,
                    'Difference (m)': difference.iloc[:10].values
                })
                print(comparison_df.to_string(index=False))
                
                # Metadata comparison
                print(f"\nMetadata comparison:")
                if meta_msl is not None and not isinstance(meta_msl, float) and station_id in meta_msl.index:
                    print(f"MSL metadata for station {station_id}:")
                    print(meta_msl.loc[station_id])
                if meta_navd is not None and not isinstance(meta_navd, float) and station_id in meta_navd.index:
                    print(f"\nNAVD metadata for station {station_id}:")
                    print(meta_navd.loc[station_id])
    
    print("\n" + "="*80)
    print("Test completed!")
    print("="*80 + "\n")

if __name__ == '__main__':
    test_datum_comparison()

