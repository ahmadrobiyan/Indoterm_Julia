"""Read raw bytes of .agg file to identify format"""
import struct

def peek(path):
    with open(path, "rb") as f:
        data = f.read(100)
        print(f"{path}:")
        print(f"  size={len(open(path,'rb').read())}")
        print(f"  hex={data[:32].hex()}")
        print(f"  ascii={[chr(b) if 32<=b<127 else '.' for b in data[:32]]}")
        # Check if it's a GEMPACK Header Array file (starts with 'GEMPAK' magic)
        if data[:6] == b'GEMPAK' or data[:6] == b'GEMPAK':
            print(f"  Format: GEMPACK Header Array")
        # Check if it's a new-format HAR (starts with BHAR or similar)
        if data[:4] in (b'BHAR', b'HAR\x00', b'\x00\x00\x00\x00'):
            print(f"  Format: New HAR format")
        print()

if __name__ == "__main__":
    base = r"C:\Users\ahmad\Downloads\INDOTERM CGE_2016_New"
    peek(f"{base}\\sec.agg")
    peek(f"{base}\\reg.agg")
    peek(f"{base}\\dashaggs\\natdis\\sec_natdis25.agg")
    peek(f"{base}\\national.har")
