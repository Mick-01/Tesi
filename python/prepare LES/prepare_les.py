import shutil
import re
import os
import stat
from typing import Optional
import pandas as pd
import numpy as np
import numpy.typing as npt

volume_path = "/Users/michele/Desktop/Tesi/foam-volume/"
run_paths = volume_path + "OpenFOAM/openfoam-v2506/run/"
datasets_paths = "/Users/michele/Desktop/Tesi/Turbulence Dataset/"
komegasst_dataset_path  = datasets_paths + "komegasst/komegasst_"
les_dataset_path        = datasets_paths + "labels/"
of_paths                = datasets_paths + "openfoam/komegasst/"

ref_df = pd.read_csv("/Users/michele/Desktop/Tesi/REF.csv")

def find_sim(sim_id:str)->str:
    match sim_id:
        case s if s.startswith("h"):
            return f"bump/{s}/"
        case s if s.startswith("case_"):
            return f"pehill/{s}/"
        case s if s.startswith("convdiv"):
            return f"convdiv/{s}/"
        case s if s.startswith("cbfs"):
            return f"cbfs/{s}/"

scenarios = {
    "BUMP": "bump/",
    "PHLL" : "pehill/",
    "CDNV" : "convdiv/convdiv",
    "CBFS" : "cbfs/cbfs",
    "DUCT" : "squareDuct/squareDuct_Re"
}

# returns path to simulation from basic string
def get_sim_addr(sim_id: str) -> str:
    scenario_id, case_id = sim_id.split("_", maxsplit=1)
    return scenarios[scenario_id] + case_id

def get_les_speed(sim_id: str) -> np.ndarray:
    df = ref_df[ref_df["Case"]==sim_id]
    return df[["REF_U_1","REF_U_2","REF_U_3"]].to_numpy()

# changes the content of the file to prepare it for frozenRANS
def change_in_file(filename: str, content:npt.ArrayLike)-> None:
    with open(filename, "r") as file:
        f = file.read()
        content_str = ""
        type_str = ""
    if len(content.shape) == 1:
        content_str = str(content.shape[0]) + "\n" \
            + "(\n"+"\n".join(content.astype(str))+"\n)"
        type_str = "nonuniform List<scalar>"
    elif len(content.shape) ==2:
        content_str = str(content.shape[0]) + "\n" \
            + "(\n"+"\n".join(["("+" ".join(c.astype(str))+")" for c in content])+"\n)"
        type_str = "nonuniform List<vector>"
    s = re.sub(r"(internalField\s+)(.*?);",f"\\1{type_str}\n{content_str};",f, re.S)
    with open(filename,"w") as file:
        file.write(s)    
    

# prepares a simulation to be run in openfoam
def prepare_simulation(sim_id: str) -> None:
    # creates paths needed afterwards
    sim_path = of_paths + find_sim(sim_id)
    run_path = run_paths + find_sim(sim_id)

    # copies everything to the run directory
    shutil.copytree(sim_path,run_path)
    # modifies to run frozenRans
    U = get_les_speed(sim_id)
    print(U.shape)
    print(len(U.T))
    U_path = run_path+"/0/U"
    k = ref_df[ref_df["Case"]==sim_id]["REF_k"]
    k_path = run_path+"/0/k"
    change_in_file(U_path,U)
    change_in_file(k_path,k)
    with open(run_path+"/constant/turbulenceProperties", "r+") as file:
        turb = file.read()
        file.seek(0)
        file.write(re.sub(r"(RASModel\s+)(.*?);",r"\1frozenSST;",turb))
    with open(run_path+"/Runsim", "r+") as file:
        script = file.read()
        file.seek(0)
        file.write(re.sub("[0-9]+ simpleFoam","4 frozenFoam",script))
    with open(run_path+"/system/decomposeParDict", "r+") as file:
        script = file.read()
        file.seek(0)
        file.write(re.sub(r"(numberOfSubdomains\s+)(.*?);",r"\1 4;",script))
    os.chmod(run_path+"/Runsim",stat.S_IRWXU | stat.S_IRWXG | stat.S_IRWXO)
    
def run_simulation(sim_id: str, dst_path: Optional[str] = None) -> None:
    cwd = os.getcwd()
    print(cwd)
    if dst_path is None:
        dst_path = cwd
    os.chdir(volume_path)
    sim_path = "OpenFOAM/openfoam-v2506/run/" + get_sim_addr(sim_id)
    command = f'''foam -c "sh -c './{sim_path}/Runsim' " '''
    os.system(command)
    # shutil.copyfile(sim_path+"/sim.log",dst_path)
    # shutil.copyfile(sim_path+"/error.log",dst_path)
    os.chdir(cwd)
    

if __name__ == "__main__":
    sim_id = "h20"
    prepare_simulation(sim_id)
    #run_simulation(sim_id)