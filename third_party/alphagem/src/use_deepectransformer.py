import os
import subprocess
import shutil
import sys


PINNED_RUNTIME_BOOTSTRAP = r"""
import hashlib
import sys

expected_sha256 = sys.argv[1]
expected_script_sha256 = sys.argv[2]
script_path = sys.argv[3]
digest = hashlib.sha256()
with open(sys.executable, "rb") as handle:
    for chunk in iter(lambda: handle.read(1024 * 1024), b""):
        digest.update(chunk)
if digest.hexdigest() != expected_sha256:
    raise RuntimeError("DeepEC Python executable SHA-256 mismatch")
if sys.version_info[:3] != (3, 7, 7):
    raise RuntimeError("DeepEC requires exact Python 3.7.7")
import torch
import transformers
if torch.__version__ != "1.7.0":
    raise RuntimeError("DeepEC requires exact torch 1.7.0")
if transformers.__version__ != "3.5.1":
    raise RuntimeError("DeepEC requires exact transformers 3.5.1")
with open(script_path, "rb") as handle:
    script_bytes = handle.read()
if hashlib.sha256(script_bytes).hexdigest() != expected_script_sha256:
    raise RuntimeError("DeepEC entrypoint SHA-256 mismatch")
sys.argv = sys.argv[3:]
namespace = {
    "__name__": "__main__",
    "__file__": script_path,
    "__package__": None,
    "__cached__": None,
}
exec(compile(script_bytes, script_path, "exec"), namespace)
""".strip()

def activate_deepectransformer_env():
    try:
        subprocess.run("source activate deepectransformer", shell=True, check=True)
        print("Switched to deepectransformer environment.")
    except subprocess.CalledProcessError as e:
        print(f"Error switching to CLEAN environment: {e}")


def use_deepectransformer(
    name,
    python_executable=None,
    expected_python_sha256=None,
    expected_entrypoint_sha256=None,
):
    path = os.getcwd()
    original_cwd = path
    tools_dir = os.path.join(path, "tools", "DeepProZyme")
    script_path = os.path.join(path, "tools", "DeepProZyme", "run_deepectransformer.py")
    input_fasta = os.path.join(path, "working", name, f"{name}_homoleft.fasta")
    output_dir = os.path.join(path, "working", name, f"{name}_deepec_result")
    result_path = os.path.join(output_dir, "DeepECv2_result.txt")

    if not os.path.isdir(tools_dir):
        raise FileNotFoundError(f"DeepProZyme directory not found: {tools_dir}")
    if not os.path.isfile(input_fasta):
        raise FileNotFoundError(f"Input FASTA not found: {input_fasta}")

    os.makedirs(output_dir, exist_ok=True)
    if expected_python_sha256 is not None and os.path.exists(result_path):
        raise FileExistsError(
            "A formal DeepEC stage will not accept a pre-existing result; "
            "resume the validated stage or use a fresh job name"
        )
    if python_executable is not None and (
        expected_python_sha256 is None
        or expected_entrypoint_sha256 is None
    ):
        raise ValueError(
            "An explicit DeepEC interpreter requires preflight SHA-256 "
            "records for both Python and the entrypoint"
        )

    def _exists(cmd: str) -> bool:
        return shutil.which(cmd) is not None

    def _candidate_cmds():
        if python_executable is not None:
            return [[str(python_executable)]]
        cmds = []
        conda_exe = os.environ.get("CONDA_EXE")
        if conda_exe and os.path.isfile(conda_exe):
            cmds.append([conda_exe, "run", "-n", "deepectransformer", "python"])
        if _exists("conda"):
            cmds.append(["conda", "run", "-n", "deepectransformer", "python"])
        if _exists("mamba"):
            cmds.append(["mamba", "run", "-n", "deepectransformer", "python"])
        if _exists("micromamba"):
            cmds.append(["micromamba", "run", "-n", "deepectransformer", "python"])
        cmds.append([sys.executable])
        return cmds

    try:
        os.chdir(tools_dir)
        #activate_deepectransformer_env()
        base_args = [
            script_path,
            "-i", input_fasta,
            "-o", output_dir,
            "-g", "cpu",
            "-b", "128",
            "-cpu", "2",
        ]
        last_err = None
        for runner in _candidate_cmds():
            cmd = (
                runner
                + [
                    "-c",
                    PINNED_RUNTIME_BOOTSTRAP,
                    expected_python_sha256,
                    expected_entrypoint_sha256,
                ]
                + base_args
                if python_executable is not None
                else runner + base_args
            )
            try:
                print("Running:", " ".join(cmd))
                subprocess.run(cmd, check=True)
                print("deepectransformer completed.")
                break
            except subprocess.CalledProcessError as e:
                last_err = e
                print(f"Runner failed: {' '.join(runner)}; return code {e.returncode}")
        else:
            if last_err:
                raise last_err
    finally:
        os.chdir(original_cwd)
    if not os.path.isfile(result_path):
        raise RuntimeError(
            f"DeepECTransformer completed without its expected result: {result_path}"
        )
