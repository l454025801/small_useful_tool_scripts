import numpy as np
import sys 

def number_to_mol(number, x, y=None, z=None):
    if (y==None and z==None):
        v = float(x)**3/10
        print("%s molecule(s) in a %s x %s x %s nm^3 box" %(number,x,x,x))
    elif (y!=None and z!=None):
        v = float(x)*float(y)*float(z)/10
        print("%s molecule(s) in a %s x %s x %s nm^3 box" %(number,x,y,z))
    else:
        print("either specify one dimension for a cubic or three dimensions for a rectangular prism")
        return None

    mol = (int(number)/6.023)/v
    print("The concentration is %.6f mol/L" %mol)
    return mol

def main():
    try:
        number_to_mol(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
    except:
        number_to_mol(sys.argv[1], sys.argv[2])

if __name__ == '__main__':
    main()
