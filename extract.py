import os
import shutil
from pathlib import Path

# Find pip cache folder automatically
cache_dir = Path(os.path.expanduser("~")) / ".cache" / "pip" / "http-v2"
if os.name == "nt":  # Windows path adjustment
    cache_dir = Path(os.getenv("LOCALAPPDATA", "")) / "pip" / "cache" / "http-v2"

output_dir = Path("./recovered_wheels")
output_dir.mkdir(exist_ok=True)

print(f"Scanning cache at: {cache_dir}")
count = 0

for file_path in cache_dir.glob("**/*"):
    if file_path.is_file() and not file_path.name.endswith(".body"):
        try:
            with open(file_path, "rb") as f:
                content = f.read()

            # Search for the zip/wheel file signature hidden inside the cache file
            wheel_start = content.find(b"PK\x03\x04")
            if wheel_start != -1:
                # Look for the filename inside the HTTP metadata
                filename = "package.whl"
                for chunk in content.split(b"\x00"):
                    if b".whl" in chunk:
                        try:
                            filename = chunk.decode("utf-8", errors="ignore").split("/")[-1]
                            if not filename.endswith(".whl"):
                                filename = filename.split("\\")[-1]
                        except:
                            pass

                # If we couldn't parse the name perfectly, give it a unique ID
                if ".whl" not in filename or len(filename) > 100:
                    filename = f"extracted_{count}.whl"

                out_file = output_dir / filename
                with open(out_file, "wb") as out_f:
                    out_f.write(content[wheel_start:])
                print(f"-> Extracted: {filename}")
                count += 1
        except Exception:
            pass

print(f"\nDone! Extracted {count} files to the './recovered_wheels' folder.")
