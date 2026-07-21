import sys
import harpy

def export(har_path, csv_path):
    hf = harpy.HarFileObj(har_path)
    names = hf.getRealHeaderArrayNames()
    with open(csv_path, "w", encoding="utf-8") as out:
        for name in names:
            obj = hf._getHeaderArrayObj(name)
            arr = obj.array
            set_names = obj.setNames
            # set element labels per dimension, in order
            elem_lists = [
                [str(e).strip() for e in dim_elems]
                for dim_elems in obj.sets.setElements
            ] if set_names else []

            ndim = arr.ndim
            if ndim == 0:
                out.write(f"__scalar__,Value\n")
                out.write(f"{name},{float(arr)}\n")
                continue

            cols = ",".join(f"dim{i+1}" for i in range(ndim))
            out.write(f"__header__:{name},{cols},Value\n")

            import numpy as np
            it = np.nditer(arr, flags=["multi_index"])
            for val in it:
                idx = it.multi_index
                labels = []
                for d, i in enumerate(idx):
                    if d < len(elem_lists) and elem_lists[d]:
                        labels.append(elem_lists[d][i])
                    else:
                        labels.append(str(i + 1))
                out.write(",".join(labels) + f",{float(val)}\n")
    print(f"Exported {len(names)} headers from {har_path} -> {csv_path}")

if __name__ == "__main__":
    har_path, csv_path = sys.argv[1], sys.argv[2]
    export(har_path, csv_path)
