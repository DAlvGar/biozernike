import subprocess

def run_command(command):
    result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    return result.returncode, result.stdout, result.stderr

def reconstruct(input_file, zernike_order, output_file, dim=32, origin=None, save_cache=False, scale=None, use_cache=True, width=None):
    args = ["ligzernike", "reconstruct", "-i", input_file, "-N", str(zernike_order), "-o", output_file, "--dim", str(dim)]
    if origin:
        args.extend(["--origin"] + [str(o) for o in origin])
    if save_cache:
        args.append("--saveCache")
    if scale:
        args.extend(["--scale", str(scale)])
    if use_cache:
        args.append("--useCache")
    if width:
        args.extend(["--width", str(width)])
    return run_command(args)

def rotate(m1, m2, zernike_order):
    return run_command(["ligzernike", "rotate", "-m1", m1, "-m2", m2, "-N", str(zernike_order)])
