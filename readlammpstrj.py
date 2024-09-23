import numpy as np
import MDAnalysis as mda
import pandas as pd
import plotly.graph_objects as go
import re

def read_lammpstrj(lammpstrj, steps, frames, trj_format):
    '''
    functions to read lammps dump trajectory. Support grand canonical ensemble (charge regulation, RXNff, MC).
    Return a dictionary containing all frames, each frame is a pandas Dataframe containing all atoms.
    
     lammpstrj: File path. lammpsdump file
         steps: Int. total steps in the simulation
        frames: Int. how many frames are saved in the trj file
    trj_format: List. the dumped file format, ex: ["index", "type", "x", "y", "z"]
    '''
    with open(lammpstrj) as lmptrj:
        content = lmptrj.read()
    rst_dict = {}
    save_frames = (np.linspace(0, steps, frames)).astype(int)
    
    # All other frames beside last one
    for frame in save_frames[:-1]:
        print("*", end='')
        re_txt = 'TIMESTEP\n%i\n.*?TIMESTEP' % frame
        text = re.search(re_txt, content, re.DOTALL).group()
        
        a = text.splitlines()[9:-1]   # split into lines, only pcik atoms and discard other lines
        b = [re.split(" ", line) for line in a] # split each line into elements
        c = np.asarray(b, dtype=float)  # convert text into numbers
        
        frame_df = pd.DataFrame(data=c, columns = trj_format)
        rst_dict[frame] = frame_df
    
    # Last frame
    re_txt = 'TIMESTEP\n%i\n.*' % save_frames[-1]
    text = re.search(re_txt, content, re.DOTALL).group()
        
    a = text.splitlines()[9:-1]   # split into lines, only pcik atoms and discard other lines
    b = [re.split(" ", line) for line in a] # split each line into elements
    c = np.asarray(b, dtype=float)  # convert text into numbers
    print(len(c))
    frame_df = pd.DataFrame(data=c, columns = trj_format)
    rst_dict[save_frames[-1]] = frame_df    
    
    lmptrj.close()
    return rst_dict

