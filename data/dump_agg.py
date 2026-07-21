"""Extract aggregation mapping from GEMPACK binary files using harpy's sl4 module"""
import harpy
import sys

def dump_agg(path, label):
    print(f"\n=== {label}: {path} ===")
    try:
        # Try sl4 (compact data) format
        sl = harpy.SL4(path)
        print(f"Type: SL4")
        print(f"Headers: {dir(sl)}")
        # The agg file might have datasets
        for name in dir(harpy.sl4):
            if not name.startswith('_'):
                print(f"  sl4.{name}")
    except:
        pass
    try:
        # Try to read as a GEMPACK data file
        hf = harpy.HarFileObj(path)
        names = hf.getRealHeaderArrayNames()
        print(f"Headers (HAR): {names}")
    except Exception as e:
        # Try direct reading
        try:
            import harpy.har_file_io as hio
            data = hio.read(path)
            print(f"har_file_io result type: {type(data)}")
            if hasattr(data, 'keys'):
                for k in data.keys():
                    print(f"  key={k}, type={type(data[k])}")
        except Exception as e2:
            print(f"  HAR error: {e}")
            print(f"  hio error: {e2}")

if __name__ == "__main__":
    base = r"C:\Users\ahmad\Downloads\INDOTERM CGE_2016_New"
    dump_agg(f"{base}\\sec.agg", "Sector aggregation (185->25)")
    dump_agg(f"{base}\\reg.agg", "Region aggregation (34->34)")
    dump_agg(f"{base}\\dashaggs\\natdis\\sec_natdis25.agg", "natdis sector agg")
    dump_agg(f"{base}\\dashaggs\\natdis\\reg_natdis34.agg", "natdis region agg")
