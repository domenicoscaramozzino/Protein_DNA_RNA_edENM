try:
    import MDAnalysis as mda
except ImportError:
    raise ImportError("Please install MDAnalysis to run this script")
import numpy as np
import sys
import scipy.linalg as linalg
import warnings
warnings.filterwarnings("ignore")

from Bio.PDB import MMCIFParser, Structure, Chain, Residue, Atom, PDBIO, Model, PDBParser
from collections import defaultdict
import os
import math


# define the base class for all Elastic Network Model classes
class BaseENM(object):
    def __init__(self, structure, selection, cutoff, gamma_function, define_spring_type):
        self._structure = structure
        self._selection = selection
        self._cutoff = cutoff
        self._gamma_function = gamma_function
        self._gamma_matrix = None
        self._define_spring_type = define_spring_type
        self._hessian = None
        self.nodeAtomNames = None
        self.nodeResNames = None
        self.nodeResIds = None
        self.nodeChainIds = None

    def _change_resname(self,resname):
        if "A" in resname and resname != "GUA" and resname != "URA":
            newname = "A"
        elif "U" in resname and resname != "GUA":
            newname = "U"
        elif "G" in resname:
            newname = "G"
        elif "C" in resname:
            newname = "C"
        elif "T" in resname:
            newname = "T"
        else:
            print(f"Error: residue {resname} not recognized")
            sys.exit(1)
        return newname

    def _get_coordinates(self):
        u = mda.Universe(self._structure).select_atoms(self._selection)
        coordinates = u.atoms.positions
        self.nodeAtomNames = np.array([atom.name for atom in u.atoms])
        self.nodeResNames = np.array([self._change_resname(atom.resname) for atom in u.atoms])
        self.nodeResIds = np.array([atom.resid for atom in u.atoms])
        self.nodeChainIds = np.array([atom.chainID for atom in u.atoms])
        return coordinates

    def build_hessian(self, sparse = False):
        coordinates = self._get_coordinates()
        n_atoms = coordinates.shape[0]
        if not sparse:
            hessian = np.zeros((n_atoms*3, n_atoms*3))
            gamma_matrix = np.zeros((n_atoms, n_atoms))
        else:
            from scipy import sparse as scipy_sparse
            hessian = scipy_sparse.lil_matrix((n_atoms*3, n_atoms*3))
            gamma_matrix = scipy_sparse.lil_matrix((n_atoms, n_atoms))
        
        
        for i in range(n_atoms):
            i3, i33 = i*3, i*3+3
            i_p1 = i+1
            i2j_all = coordinates[i_p1:] - coordinates[i]
            for j, dist2 in enumerate(np.sum(i2j_all**2, axis=1)):
                if dist2 <= self._cutoff**2:
                    i2j = i2j_all[j]
                    j += i_p1
                    j3, j33 = j*3, j*3+3
                    bond_type = self._define_spring_type(self.nodeAtomNames[i], self.nodeResIds[i], self.nodeResNames[i],
                                            self.nodeAtomNames[j], self.nodeResIds[j], self.nodeResNames[j],
                                            self.nodeChainIds[i], self.nodeChainIds[j], np.sqrt(dist2))
                    gamma = self._gamma_function(bond_type, np.sqrt(dist2))
                    super_element = np.outer(i2j, i2j) * (-gamma/dist2)
                    hessian[i3:i33, j3:j33] = super_element
                    hessian[j3:j33, i3:i33] = super_element
                    hessian[i3:i33, i3:i33] -= super_element
                    hessian[j3:j33, j3:j33] -= super_element
                    gamma_matrix[i,j] = gamma
                    gamma_matrix[j,i] = gamma
        if sparse:
            hessian = hessian.tocsr()
        self._hessian = hessian
        self._gamma_matrix = gamma_matrix

    def calcModes(self, keepzeros=False, sparse=False, k=None):
        if sparse:
            if k is None:
                k = self._hessian.shape[0]-1
            elif not k < self._hessian.shape[0]:
                raise ValueError('Using the sparse option it is possible to compute a number of modes '
                                 'lower than the degrees of freedom of the system.')
            from scipy.sparse import issparse
            if not issparse(self._hessian):
                eigvals, eigvecs = linalg.eigh(self._hessian)
            else:
                from scipy.sparse import linalg as scipy_sparse_la
                eigvals, eigvecs = scipy_sparse_la.eigsh(self._hessian, k=k, which='SA')
        else:
            eigvals, eigvecs = linalg.eigh(self._hessian)
        n_zeros = np.sum(eigvals < 1e-10)
        if n_zeros > 6:
            print(f"Warning: {n_zeros} zero eigenvalues found. Check the structure and the selection.")
        if not keepzeros:
            eigvals = eigvals[6:]
            eigvecs = eigvecs[:,6:]
        return eigvals, eigvecs    

class BussiENM(BaseENM):
    """
    class to implement the nucleic acid ENM from Giovanni Pinamonti, Sandro Bottaro, Cristian Micheletti and Giovanni Bussi
    from the paper:
    Elastic network models for RNA: a comparative assessment with molecular dynamics and SHAPE experiments
    doi: 10.1093/nar/gkv708

    Features of the ENM
    nodes: atom names P C1' and C2
    cutoff: 9 Angstroms
    gamma: 1.0 between every pair of nodes
    """
    def __init__(self, structure, selection="name P C1' C2", cutoff=9.0):
        super().__init__(structure, selection, cutoff, self._gamma, self._spring_type)

    def _gamma(self, *args):
        # function to return the gamma value for the spring constant
        return 1.0
    
    def _spring_type(self, *args):
        # dummy function to return the bond type introduced for compatibility with the parent class
        return 0
        
class NucleicAcidENM(BaseENM):
    """
    New class to implement the nucleic acid ENM developed by Marco Cannariato, Domenico Scaramozzino, Byung Ho Lee, and Laura Orellana
    during the exchange period of Marco Cannariato at the Protein Dynamics and Mutation Lab at the Karolisnka Institute in Stockholm, Sweden
    """
    def __init__(self, structure, selection="name P C1' C2"):
        super().__init__(structure, selection, self._define_cutoff(structure), self._get_gamma, self._spring_type) # 

    def _define_cutoff(self, structure):
        return 11 
    
    def _spring_type(self,atom1, resid1, resname1, atom2, resid2, resname2, chain1, chain2, r):
        # format the residue names in simple way
        couple_res = f"{resname1}-{resname2}"
        couple_atom = f"{atom1}-{atom2}"
        if (resid1 == resid2) & (couple_atom in ["C1'-C2","C2-C1'"]):
                return "covalent_pir" if resname1 in ["C","U","T"] else "covalent_pur"
        elif (couple_res in ["C-G","G-C"]) & (couple_atom == "C2-C2") & (r >= 4.1) & (r <= 4.4):
            return "hbond"
        elif (couple_res == "G-G") & (couple_atom == "C2-C2") & (r >= 6.7) & (r <= 7.0):
            return "hbond"
        elif (couple_atom in ["P-P","P-C1'","C1'-P"]) & (chain1 == chain2):
            return "pseudocov"
        else:
            return "vdw"

    def _get_gamma(self,bond_type,r):
        if bond_type == "hbond":
            return 65
        elif bond_type == "vdw":
            return (25/r)**1
        elif bond_type == "pseudocov":
            return (20/r)**2.8
        elif bond_type == "covalent_pir":
            return 290
        elif bond_type == "covalent_pur":
            return 120
        else:
            print(f"Error: bond type {bond_type} not recognized")
            sys.exit(1)

def calcVariances(eigvals):
    return 1/eigvals

def calcExplainedVariance(eigvals):
    return (1/eigvals)/np.sum(1/eigvals)

def calcFluctuations(eigvecs, eigvals, n_modes):
    eigenvectors = eigvecs[:,:n_modes].reshape(-1, 3, n_modes)
    variances = 1/eigvals[:n_modes]
    return np.sqrt((linalg.norm(eigenvectors,axis=1)**2/variances).sum(axis=1))

def calcOverlap(enm_vecs, pca_vecs, n_enm, n_pca=None):
    if n_pca is None:
        n_pca = n_enm
    # ensure that the vectors are normalized
    if not np.allclose(np.linalg.norm(enm_vecs, axis=0), 1):
        enm_vecs = enm_vecs/np.linalg.norm(enm_vecs, axis=0)
    if not np.allclose(np.linalg.norm(pca_vecs, axis=0), 1):
        pca_vecs = pca_vecs/np.linalg.norm(pca_vecs, axis=0)
    overlaps = np.abs(np.matmul(enm_vecs[:,:n_enm].T,pca_vecs[:,:n_pca]))
    if n_pca == 1:
        overlaps = overlaps.squeeze()
    cum_overlaps = np.sum(overlaps**2, axis=0)
    return overlaps, cum_overlaps

def calcRMSIP(enm_vecs, pca_vecs, n_modes):
    overlaps, _ = calcOverlap(enm_vecs, pca_vecs, n_modes)
    return np.sqrt(np.sum(overlaps**2)/n_modes)

def calcCollectivity(eigvecs):
    # reshape the eigenvectors to a matrix with shape (n_atoms, 3, n_modes)
    u2 = np.sum(eigvecs.reshape(-1, 3 ,eigvecs.shape[1])**2, axis=1)
    alpha = 1/np.sum(u2,axis=0)
    return np.exp(-np.sum((alpha*u2)*np.log(alpha*u2),axis=0))/u2.shape[0]

def makeTraj(eigvecs, structure, output="ENM_trj.pdb",
             steps=5, ampl=1, modes=1, selection = "name P"):
    """
    Function to create a trajectory of motion starting from a set of eigenvectors.
    """
    if isinstance(eigvecs,list):
        eigvecs = np.column_stack(eigvecs)
    if isinstance(modes, int):
        modes = [modes]
    n_modes = len(modes)
    if n_modes > eigvecs.shape[1]:
        raise ValueError("Number of modes requested is higher than the number of modes available")
    init_pos = mda.Universe(structure).select_atoms(selection).positions
    
    #n_frames=steps*2
    #xx = np.arange(-1*steps,steps+1)
    #angles = np.pi * (xx[xx!=0] + steps) / n_frames
    
    n_frames = steps * 4
    angles = np.linspace(0, 2 * np.pi, n_frames)
    
    # Use the sine of the angle to modulate the factor
    weights = np.sin(angles) * ampl
    # weights = np.linspace(-1,1,n_frames)
    for i in range(n_modes):
        projected = init_pos.reshape(-1,1) + np.outer(weights,eigvecs[:,i]).T
        new_coords = projected.T.reshape(n_frames,-1,3)
        u = mda.Merge(mda.Universe(structure).select_atoms(selection)).load_new(new_coords, order="fac")
        filename = output.split(".")[0] + f"_mode-{i+1}.pdb"
        with mda.Writer(filename, multiframe=True) as W:
            for _ in u.trajectory:
                W.write(u)

class NormalModeAnalysis2:
    def __init__(self, work_dir=None, script_dir = None,
                 cutoff_inter = 10,
                 p_high = 1.3,
                 C_low = 35,
                 C_high = 10,
                 gamma = 1,
                 dr = 8,
                 power = 0,
                 modes=10):
       
        self.work_dir = work_dir
        self.script_dir = script_dir
        self.gamma = float(gamma)
        self.dr = float(dr)
        self.power = float(power)
        self.hessian_matrices = {}
        self.coord_dict = {}
        self.eigenvalues = {}
        self.eigenvectors = {}
        self.modes = modes
        self.masses = {}
        self.hessian_weighted_matrices = {}
        self.new_positions_after_NMA = {} # contains the endpoints of the vectors that are used to map the normal modes 
        self.eigenvecs_reshaped = {}
        self.chain_entity_map = {}
        self.gamma_matrices = {}
        self.cutoff_inter = cutoff_inter
        self.p_high = p_high
        self.C_low = C_low
        self.C_high = C_high

    def _change_resname(self,resname):
        if "A" in resname and resname != "GUA" and resname != "URA":
            newname = "A"
        elif "U" in resname and resname != "GUA":
            newname = "U"
            
        elif "G" in resname:
            newname = "G"
        elif "C" in resname:
            newname = "C"
        elif "T" in resname:
            newname = "T"
        else:
            print(f"Error: residue {resname} not recognized")
            sys.exit(1)
        return newname
   
    def extract_ca_coordinates_and_masses(self, pdb, name=None):
        #helper function to parse a pdb file and  

        lst =  [('VAL',99), ('ILE',113), ('LEU',113), ('GLU',129), ('GLN',128),
                        ('ASP',115), ('ASN',114), ('HIS',137), ('TRP',186), ('PHE',147), ('TYR',163), 
                        ('ARG',156), ('LYS',128), ('SER',87), ('THR',101), ('MET',131), ('ALA',71), 
                        ('GLY',57), ('PRO',97), ('CYS',103)]
        
        aa_dict = defaultdict(lambda: 100, lst) # 100 if non canonical
        protein_residues = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
                            "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]
        nucleic_residues = ["DA", "DC", "DG", "DT", "A", "C", "G", "U", "T", "RA", "RC", "RG",
                            "RU", "ADE", "CYT", "GUA", "THY", "URI"]

        if pdb.endswith(".pdb"):
            parser = PDBParser()
        else:
            parser = MMCIFParser()
        
        if name is None:
            name = pdb.split("/")[-1].split(".")[0]
        
        structure = parser.get_structure("protein", pdb) # modify name ?
        
        coords = []
        chain_lbl = []
        masses = []
        resnames = []
        atomnames = []
        resnums = []
        
        for model in structure:
            for chain in model:
                # check if chain is protein
                if np.any([residue.get_resname() in protein_residues for residue in chain]):
                    self.chain_entity_map[chain.get_id()] = "protein"
                    for residue in chain:
                        for atom in residue:
                            if atom.get_name() == "CA":
                                atom_coords = atom.get_coord()
                                coords.append(atom_coords) # coord
                                chain_lbl.append(chain.get_id()) # chain
                                res = residue.get_resname()
                                resnames.append(res) # save the resname for protein-nucleic interaction
                                
                                res_weight = np.float64(aa_dict[res]) #* self.mass_conv_u_to_kg # convert to kg from Dalton/Unit

                                masses.append(res_weight) # mass_in_kg
                                atomnames.append("CA")
                                resnums.append(residue.get_id()[1])
                elif np.any([residue.get_resname() in nucleic_residues for residue in chain]):
                    self.chain_entity_map[chain.get_id()] = "nucleic"
                    for residue in chain:
                        tmp = {"P": [], "C1'": [], "C2": []} # in this way we are sure to have P-C1-C2 order -> needed for the ENM
                        for atom in residue:
                            if atom.get_name() in tmp.keys():
                                atom_coords = atom.get_coord()
                                tmp[atom.get_name()] = [atom_coords,
                                                        chain.get_id(),
                                                        self._change_resname(residue.get_resname()), # save resname for definition of the spring constant
                                                        np.float64(aa_dict[residue.get_resname()])/3,
                                                        atom.get_name(),
                                                        residue.get_id()[1]]
                        for v in tmp.values():
                            if v:
                                coords.append(v[0])
                                chain_lbl.append(v[1])
                                resnames.append(v[2])
                                masses.append(v[3])
                                atomnames.append(v[4])
                                resnums.append(v[5])
                else:
                    continue

        self.coord_dict[name] = (np.array(coords), chain_lbl, resnames, atomnames, resnums)
        self.masses[name] = masses

    def _find_protein_nucleic_contacts(self, pdb_code,cutoff):
        """
        Function to find protein-nucleic contacts in a structure. In this model, only one spring is
        set between each protein residue and nucleic acid residue, which is between the C-alpha atom
        of the protein residue and the closest bead of the nucleic acid residue.
        The function returns the indices of the couple of interacting residues, shape (n_contacts, 2).
        """
        # Extract the coordinates of the C-alpha atoms of the protein and the P, C1', and C2 atoms of the nucleic acids
        coords, chain_lbl, _, atomnames, resnums = self.coord_dict[pdb_code]
        protein_coords = coords[np.array(atomnames) == "CA"]
        protein_idx = np.where(np.array(atomnames) == "CA")[0]
        nucleic_coords = coords[np.isin(atomnames, ["P", "C1'", "C2"])]
        nucleic_idx = np.where(np.isin(atomnames, ["P", "C1'", "C2"]))[0]
        # Initialize the array to store the indices of the interacting residues
        contacts = []
        # Loop over the protein residues
        for i, protein_coord in enumerate(protein_coords):
            # Compute the distances between the protein residue and all nucleic acid residues
            distances = np.linalg.norm(nucleic_coords - protein_coord, axis=1)
            # if the lowest distance is more than the cutoff, skip
            if np.min(distances) > cutoff:
                continue
            selected_idx = []
            selected_res = []
            for k in np.argsort(distances):
                if distances[k] > cutoff:
                    break
                if (chain_lbl[nucleic_idx[k]],resnums[nucleic_idx[k]]) not in selected_res:
                    selected_idx.append(k)
                    selected_res.append((chain_lbl[nucleic_idx[k]],resnums[nucleic_idx[k]]))
            for k in selected_idx:
                contacts.append([protein_idx[i], nucleic_idx[k]])



            # Find the bead with the minimum distance for each
            # min_distance_idx = np.argmin(distances)
            # Store the indices of the interacting residues
            # contacts.append([protein_idx[i], nucleic_idx[min_distance_idx]])
        return contacts
       
    def construct_hessian(self, pdb_code, build_gamma=False, verbose=False):
        
        if self.coord_dict[pdb_code]:
            coords, chain_lbl, resnames, atomnames, _ = self.coord_dict[pdb_code]

        # print(coords.shape)
        n_atoms = len(coords)
        # print(n_atoms)
        hessian = np.zeros((n_atoms*3, n_atoms*3), float)
        if build_gamma:
            gamma_matrix = np.zeros((n_atoms, n_atoms), float)
        # distance_mat = np.ones((n_atoms*3, n_atoms*3), float)

        #cutoff
        # edENM_cutoff = int(2.9*np.log(n_atoms) - 2.9) # as in Lauras and Domenicos code.
        
        # if edENM_cutoff >= 20:
        #     edENM_cutoff = 20
            
        # if edENM_cutoff <= 8:
        #     edENM_cutoff = 8
        n_protein_atoms = len([c for c in chain_lbl if self.chain_entity_map[c] == "protein"])
        edENM_cutoff = int(np.clip(2.9*np.log(n_protein_atoms) - 2.9, 8, 20))
        edENM_seq_exp = 2 # so we treat neighbours different in sequence.

        n_nucleic_atoms = len([c for c in chain_lbl if self.chain_entity_map[c] == "nucleic"])
        edENM_cutoff_nucleic = 11

        edENM_cutoff_inter = self.cutoff_inter
        if verbose:
            print("Number of protein atoms: ", n_protein_atoms)
            print("Number of nucleic atoms: ", n_nucleic_atoms)

        edENM_cart_exp = 6 
        edENM_seq_force_constant = 60
        edENM_cart_force_constant = 6
        edENM_M = 3
        
        if n_nucleic_atoms > 0:
            inter_contacts = self._find_protein_nucleic_contacts(pdb_code,edENM_cutoff_inter)
            if verbose:
                print("Number of protein-nucleic spring: ", len(inter_contacts))
            if len(inter_contacts) < 1:
                print("Error! No protein-nucleic contacts found. Check the structure.")
                sys.exit(1)

        for i in range(len(coords)):
            diff = coords[i+1:, :] - coords[i] # diff between x y z of entry n and entry n-1
            squared_diff = diff**2 # squared diff for each x y and z
            # diff = diff between current x y z and  all other x y z later in sequence. 
            
            for j, s_ij in enumerate(squared_diff.sum(1)): # col wise summation   
                # j is the position of the paired bead relative to i (sequential distance is j+1)
                # s_ij is the squared distance between the two beads
                
                # check if protein-protein, protein-nucleic, or nucleic-nucleic and define the spring constant accordingly
                entity_couple = f"{self.chain_entity_map[chain_lbl[i]]}-{self.chain_entity_map[chain_lbl[j+i+1]]}"
                if entity_couple == "protein-protein":
                    # this is the normal edENM
                    if j + 1 <= edENM_M: 
                        gamma = edENM_seq_force_constant/(j+1)**edENM_seq_exp # then we will use this spring constant.
                    elif s_ij <= edENM_cutoff**2:
                        gamma = (edENM_cart_force_constant/np.sqrt(s_ij))**edENM_cart_exp 
                    else:
                        # we are ouside the cutoff, gamma = 0 so we can skip because the hessian is already initialized to 0
                        continue
                elif entity_couple == "nucleic-nucleic":
                    if s_ij <= edENM_cutoff_nucleic**2:
                        couple_res = f"{resnames[i]}-{resnames[j+i+1]}"
                        couple_atom = f"{atomnames[i]}-{atomnames[j+i+1]}"
                        if (couple_atom == "C1'-C2") & (j == 0): # covalent bond
                            gamma = 290 if resnames[i] in ["C","U","T"] else 120 # C,U,T are pirimidines
                        elif (couple_res in ["C-G","G-C"]) & (couple_atom == "C2-C2") & (np.sqrt(s_ij) > 4.1) & (np.sqrt(s_ij) < 4.4): # hbond
                            gamma = 65
                        elif (couple_res == "G-G") & (couple_atom == "C2-C2") & (np.sqrt(s_ij) > 6.7) & (np.sqrt(s_ij) < 7.0): # hbond
                            gamma = 65
                        elif (couple_atom in ["P-P","P-C1'","C1'-P"]) & (chain_lbl[i]==chain_lbl[j+i+1]): # pseudo covalent bond
                            gamma = (20/np.sqrt(s_ij))**2.8
                        else: # vdw
                            gamma = 25/np.sqrt(s_ij)
                    else:
                        # we are ouside the cutoff, gamma = 0 so we can skip because the hessian is already initialized to 0
                        continue
                elif s_ij <= edENM_cutoff_inter**2:
                    protein_idx = i if self.chain_entity_map[chain_lbl[i]] == "protein" else j+i+1
                    nucl_idx = i if self.chain_entity_map[chain_lbl[i]] == "nucleic" else j+i+1
                    if [protein_idx, nucl_idx] not in inter_contacts:
                        continue
                    # protein-nucleic
                    nucl_bead_name = atomnames[i] if self.chain_entity_map[chain_lbl[i]] == "nucleic" else atomnames[j+i+1]
                    if nucl_bead_name == "C1'":
                        gamma = self.C_low/np.sqrt(s_ij) #35
                    else:
                        protein_resname = resnames[i] if self.chain_entity_map[chain_lbl[i]] == "protein" else resnames[j+i+1]
                        protein_restype = "CN" if protein_resname in ["ASP","GLU"] else "CP" if protein_resname in ["LYS","ARG","HIS","HIP",'HIE',"HID"] \
                                else "P" if protein_resname in ["SER","THR","ASN","GLN","CYS"] else "H"
                        if ((nucl_bead_name=="C2") & (protein_restype=="CN")) | ((nucl_bead_name=="P") & (protein_restype in ["P","CP"])):
                            gamma = (self.C_high/np.sqrt(s_ij))**self.p_high #1.3
                        else:
                            gamma = self.C_low/np.sqrt(s_ij) #35
                else:
                    continue
                
                # once the spring constant is defined, we can compute the derivative and update the hessian matrix
                diff_coords = diff[j]
                j = j + i + 1
                derivative = np.outer(diff_coords, diff_coords)*(float(-gamma)/np.sqrt(s_ij)**(2+self.power))  # delta coords * -gamma/dist**0 in our case  in the paper its -gamma * delta diff  / dist**2
                hessian[i*3:i*3+3, j*3:j*3+3] = derivative
                hessian[j*3:j*3+3, i*3:i*3+3] = derivative #symmetry
                hessian[i*3:i*3+3, i*3:i*3+3] = hessian[i*3:i*3+3, i*3:i*3+3] - derivative
                hessian[j*3:j*3+3, j*3:j*3+3] = hessian[j*3:j*3+3, j*3:j*3+3] - derivative
                if build_gamma:
                    gamma_matrix[i, j] = gamma
                    gamma_matrix[j, i] = gamma
        self.hessian_matrices[pdb_code] = hessian
        if build_gamma:
            self.gamma_matrices[pdb_code] = gamma_matrix
        
    

    def add_masses_to_hessian(self, pdb_code):

        if self.hessian_matrices[pdb_code].size > 0: # not empty
            H = self.hessian_matrices[pdb_code]
            masses = self.masses[pdb_code] # masses per residue.
            N_x3 = H.shape[0] # is symmetrical and 3x num_atoms


        MWH = np.zeros_like(H, dtype=np.float64)
        mass_res_expanded = np.repeat(masses, 3) # same mass for x, y , z
        for i in range(N_x3):
            for j in range(N_x3):
                MWH[i, j] = H[i,j] / np.sqrt(mass_res_expanded[i//3] * mass_res_expanded[j//3])
                
            
        self.hessian_weighted_matrices[pdb_code] = MWH
    
    def solve_eigenproblem(self, pdb_code, weighted=False):
        # given that H is symmetrical, we can leverage eigh from numpy.linalg.

        #check symmetry
        if not self.hessian_matrices[pdb_code].size > 0:
            print("Hessian matrix is needed")
            return
    
        # Symmetrize the Hessian matrix
        if weighted:
            if self.hessian_weighted_matrices[pdb_code].size > 0:
                hessian = self.hessian_weighted_matrices[pdb_code]
        else:
            hessian = self.hessian_matrices[pdb_code]
        
        #print(f"{hessian=}")
        #with open("/home/micnag/bioinformatics/domenico_nma/test_weighted_hessian.txt", "w") as fh_out:
        #    np.savetxt(fh_out, hessian, fmt='%7.4f')
        
        #hessian_sym = (hessian + hessian.T) / 2
    
        # Solve the eigenproblem
        eigenvalues, eigenvectors = np.linalg.eigh(hessian)

        num_non_interesting_modes = 6  # Assuming the first 6 modes are not interesting
        eigenvalues = eigenvalues[num_non_interesting_modes:num_non_interesting_modes+self.modes]
        eigenvectors = eigenvectors[:, num_non_interesting_modes:num_non_interesting_modes+self.modes]
    
        # Store the results
        self.eigenvalues[pdb_code] = eigenvalues
        self.eigenvectors[pdb_code] = eigenvectors
        # now lets apply mass weighting as well

    def compute_vectors_for_NMA(self, pdb_code, num_modes=10):
        # we need to compute the vectors for each residue in the structure based on x y z triple pairs in n modes:
        # e.g if dim = 994*3, 10 we need to split first dim into pieces of 3, compute the x y z direction.
        # then we need to grab the corresponding CA residue. Map this vector ontop of the residue. 
        # and we end up with another point that is the endpoint of the vector.

        if self.eigenvectors and self.eigenvalues:
            
            self.new_positions_after_NMA[pdb_code] = {} # empty dict . for each normal mode we store results.

            CA_coords, _, _, _, _ = self.coord_dict[pdb_code] # coords and chain for the struc.


            tmp_dict = {}
            for mode in range(num_modes):
                eigenvecs_mode = self.eigenvectors[pdb_code][:, mode] # gives back the array of shape n_residues*3 , eigenvalues
                
                vecs_reshaped = eigenvecs_mode.reshape((-1, 3)) # shape (N,3) 
                
                new_positions = CA_coords + vecs_reshaped

                tmp_dict[mode] = vecs_reshaped                
                self.new_positions_after_NMA[pdb_code][f"mode_{mode+1}"] = (CA_coords, new_positions)

            self.eigenvecs_reshaped[pdb_code] = tmp_dict
        else:
            # not present eigvecs and eigvals
            print("Eigenvectors or eigenvalues not found for the PDB code provided.")
    
    def write_pymol_script_for_visualization(self, pdb_code, num_modes=10, scale_factor=75):

        # Define parameters for the cone and cylinder (arrow representation)
        cone_radius = 0.5  # Radius of the base of the cone, define this before using it
        cone_length = 0.8  # Length of the cone
        cylinder_radius = 0.2  # Radius of the cylinder, increased for a thicker appearance

        #print(self.new_positions_after_NMA)
        # might be problematic if num modes asked for is < then present in self.new_positions_after_NMA. needsa  check 
        script_location = f"{self.work_dir}/{os.path.basename(pdb_code[:-4])}_normal_mode_{1}_{num_modes}.pml"

        if os.path.exists(script_location):
            os.remove(script_location)

        for mode in range(num_modes):
            
            CA_coords, new_positions = self.new_positions_after_NMA[pdb_code][f"mode_{(mode+1)}"]
            
            with open(script_location, 'a') as script_file:
                script_file.write(f"load {pdb_code}, {os.path.basename(pdb_code)}\n")
                script_file.write(f"hide everything,  {os.path.basename(pdb_code)}\n")
                script_file.write(f"show cartoon,  {os.path.basename(pdb_code)}\n")
        
                for i, (ca, endpoint) in enumerate(zip(CA_coords, new_positions), start=1):
                    scaled_vector = [(e - c) * scale_factor for e, c in zip(endpoint, ca)]
                    scaled_endpoint = [c + v for c, v in zip(ca, scaled_vector)]
        
                    # Begin Python block
                    script_file.write("python\n")
                    script_file.write("from pymol.cgo import CYLINDER, CONE  # Import CGO primitives\n")
                    script_file.write("from pymol import cmd\n")
                    script_file.write("\n")
                    script_file.write("cgo = [\n")
        
                    # Define the cylinder (arrow tail) in yellow
                    script_file.write(f"    CYLINDER, {ca[0]}, {ca[1]}, {ca[2]}, {scaled_endpoint[0]}, {scaled_endpoint[1]}, {scaled_endpoint[2]}, {cylinder_radius}, 1.0, 1.0, 0.0, 1.0, 1.0, 0.0,\n")
        
                    # Define the full cone (arrow head) in yellow with a sharp tip (radius 0.0)
                    cone_tip = [scaled_endpoint[0] + scaled_vector[0] * cone_length,
                                scaled_endpoint[1] + scaled_vector[1] * cone_length,
                                scaled_endpoint[2] + scaled_vector[2] * cone_length]
                    script_file.write(f"    CONE, {scaled_endpoint[0]}, {scaled_endpoint[1]}, {scaled_endpoint[2]}, {cone_tip[0]}, {cone_tip[1]}, {cone_tip[2]}, {cone_radius}, 0.0, 1.0, 1.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0,\n")
        
                    # End CGO definition
                    script_file.write("]\n")
        
                    # Load the CGO object
                    script_file.write(f"cmd.load_cgo(cgo, 'vector_{i}')\n")
                    script_file.write("python end\n\n")
        
                script_file.write(f"zoom {os.path.basename(pdb_code)}\n")
    
        print(f"PyMOL script generated at: {script_location}")


    def animate_NMA(self, template=None, steps=5, both_sides=True, A=150, verbose=False, outpath='.'):
        # go 5 steps in both directions from eq. position.. multiply 1.1 1.2 1.3 1.4 1.5 and minus towards the direction of NM. 
        # overall quality CSO.
        # specific mode.  go 1vs1 inner mode. # abs for range [-1,1] so abs will do scaling between [0,1]

        for mmcif_code in self.coord_dict:
            # coords, chains, resnames, atomnames = self.coord_dict[mmcif_code]  # Extract original coordinates and chain information
            modes_dict = self.eigenvecs_reshaped[mmcif_code]  # Extract eigenvectors for each mode
            if template is None:
                template = mmcif_code
            if template.endswith(".pdb"):
                parser = PDBParser()
            else:
                parser = MMCIFParser()
            original_structure = parser.get_structure("template", template)

            for mode, mode_eigenvecs in modes_dict.items():
                # get only indices for C-alpha atoms and P atoms
                mode_eigenvecs = mode_eigenvecs[(np.array(self.coord_dict[mmcif_code][3]) == "P")|(np.array(self.coord_dict[mmcif_code][3]) == "CA")]

                animated_structure = Structure.Structure(f'animated_mode_{mode}')

                total_steps = steps * 2 if both_sides else steps  # Total steps considering both directions if both_sides is True
                for step in range(-1*steps, steps + 1):
                    if step == 0:
                        continue  # Skip the equilibrium position to avoid duplication

                    ## Calculate the angle for the sine function
                    #angle = math.pi * step / total_steps
                    ## Use the sine of the angle to modulate the factor
                    #factor = math.sin(angle) * A

                    angle = 2 * math.pi * (step / total_steps)   # ONE FULL sine oscillation
                    factor = math.sin(angle) * A

                    new_model = Model.Model(step + steps)  # Adjust model ID to be positive

                    ca_atom_index = 0  # Initialize C-alpha atom index
                    for orig_chain in original_structure.get_chains():
                        new_chain = Chain.Chain(orig_chain.id)

                        for orig_residue in orig_chain.get_residues():
                            for orig_atom in orig_residue:
                                if orig_atom.get_name() in ['CA','P']:  # Filter for C-alpha atoms
                                    if ca_atom_index >= len(mode_eigenvecs):  # Check for index out of bounds
                                        break
                                    eigenvec = mode_eigenvecs[ca_atom_index]
                                    original_coord = np.array(orig_atom.coord)
                                    delta = eigenvec * factor
                                    new_coord = original_coord + delta

                                    new_atom = Atom.Atom(orig_atom.name, new_coord, orig_atom.bfactor, orig_atom.occupancy, orig_atom.altloc, orig_atom.fullname, orig_atom.serial_number, orig_atom.element)
                                    new_residue = Residue.Residue(orig_residue.id, orig_residue.resname, orig_residue.segid)
                                    new_residue.add(new_atom)
                                    new_chain.add(new_residue)

                                    ca_atom_index += 1  # Increment C-alpha atom index

                        new_model.add(new_chain)
                    animated_structure.add(new_model)

                io = PDBIO()
                io.set_structure(animated_structure)
                if (mmcif_code.endswith(".pdb")) | (mmcif_code.endswith(".cif")):
                    output_file_path = f"{outpath}/{mmcif_code[:-4]}_mode_{mode+1}_animated.pdb"
                else:
                    output_file_path = f"{outpath}/{mmcif_code}_mode_{mode+1}_animated.pdb"
                io.save(output_file_path)
                if verbose:
                    print(f"Animated PDB file for mode {mode+1} saved as {output_file_path}")



class edENM:
    def __init__(self, work_dir=None, script_dir = None,
                 cutoff_inter = 10,
                 p_high = 1.3,
                 C_low = 35,
                 C_high = 10,
                 gamma = 1,
                 dr = 8,
                 power = 0,
                 modes=10):
       
        self.work_dir = work_dir
        self.script_dir = script_dir
        self.gamma = float(gamma)
        self.dr = float(dr)
        self.power = float(power)
        self.hessian_matrices = {}
        self.coord_dict = {}
        self.eigenvalues = {}
        self.eigenvectors = {}
        self.modes = modes
        self.masses = {}
        self.hessian_weighted_matrices = {}
        self.new_positions_after_NMA = {} # contains the endpoints of the vectors that are used to map the normal modes 
        self.eigenvecs_reshaped = {}
        self.chain_entity_map = {}
        self.gamma_matrices = {}
        self.cutoff_inter = cutoff_inter
        self.p_high = p_high
        self.C_low = C_low
        self.C_high = C_high

    def _change_resname(self,resname):
        if "A" in resname and resname != "GUA" and resname != "URA":
            newname = "A"
        elif "U" in resname and resname != "GUA":
            newname = "U"
        elif "G" in resname:
            newname = "G"
        elif "C" in resname:
            newname = "C"
        elif "T" in resname:
            newname = "T"
        else:
            print(f"Error: residue {resname} not recognized")
            sys.exit(1)
        return newname
   
    def extract_ca_coordinates_and_masses(self, pdb, name=None):
        #helper function to parse a pdb file and  

        lst =  [('VAL',99), ('ILE',113), ('LEU',113), ('GLU',129), ('GLN',128),
                        ('ASP',115), ('ASN',114), ('HIS',137), ('TRP',186), ('PHE',147), ('TYR',163), 
                        ('ARG',156), ('LYS',128), ('SER',87), ('THR',101), ('MET',131), ('ALA',71), 
                        ('GLY',57), ('PRO',97), ('CYS',103)]
        
        aa_dict = defaultdict(lambda: 100, lst) # 100 if non canonical
        protein_residues = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
                            "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]
        nucleic_residues = ["DA", "DC", "DG", "DT", "A", "C", "G", "U", "T", "RA", "RC", "RG",
                            "RU", "ADE", "CYT", "GUA", "THY", "URI"]

        if pdb.endswith(".pdb"):
            parser = PDBParser()
        else:
            parser = MMCIFParser()
        
        if name is None:
            name = pdb.split("/")[-1].split(".")[0]
        
        structure = parser.get_structure("protein", pdb) # modify name ?
        
        coords = []
        chain_lbl = []
        masses = []
        resnames = []
        atomnames = []
        resnums = []
        
        for model in structure:
            for chain in model:
                # check if chain is protein
                if np.any([residue.get_resname() in protein_residues for residue in chain]):
                    self.chain_entity_map[chain.get_id()] = "protein"
                    for residue in chain:
                        for atom in residue:
                            if atom.get_name() == "CA":
                                atom_coords = atom.get_coord()
                                coords.append(atom_coords) # coord
                                chain_lbl.append(chain.get_id()) # chain
                                res = residue.get_resname()
                                resnames.append(res) # save the resname for protein-nucleic interaction
                                
                                res_weight = np.float64(aa_dict[res]) #* self.mass_conv_u_to_kg # convert to kg from Dalton/Unit

                                masses.append(res_weight) # mass_in_kg
                                atomnames.append("CA")
                                resnums.append(residue.get_id()[1])
                elif np.any([residue.get_resname() in nucleic_residues for residue in chain]):
                    self.chain_entity_map[chain.get_id()] = "nucleic"
                    for residue in chain:
                        tmp = {"P": [], "C1'": [], "C2": []} # in this way we are sure to have P-C1-C2 order -> needed for the ENM
                        for atom in residue:
                            if atom.get_name() in tmp.keys():
                                atom_coords = atom.get_coord()
                                tmp[atom.get_name()] = [atom_coords,
                                                        chain.get_id(),
                                                        self._change_resname(residue.get_resname()), # save resname for definition of the spring constant
                                                        np.float64(aa_dict[residue.get_resname()])/3,
                                                        atom.get_name(),
                                                        residue.get_id()[1]]
                        for v in tmp.values():
                            if v:
                                coords.append(v[0])
                                chain_lbl.append(v[1])
                                resnames.append(v[2])
                                masses.append(v[3])
                                atomnames.append(v[4])
                                resnums.append(v[5])
                else:
                    continue

        self.coord_dict[name] = (np.array(coords), chain_lbl, resnames, atomnames, resnums)
        self.masses[name] = masses

    def _find_protein_nucleic_contacts(self, pdb_code,cutoff):
        """
        Function to find protein-nucleic contacts in a structure. In this model, only one spring is
        set between each protein residue and nucleic acid residue, which is between the C-alpha atom
        of the protein residue and the closest bead of the nucleic acid residue.
        The function returns the indices of the couple of interacting residues, shape (n_contacts, 2).
        """
        # Extract the coordinates of the C-alpha atoms of the protein and the P, C1', and C2 atoms of the nucleic acids
        coords, chain_lbl, _, atomnames, resnums = self.coord_dict[pdb_code]
        protein_coords = coords[np.array(atomnames) == "CA"]
        protein_idx = np.where(np.array(atomnames) == "CA")[0]
        nucleic_coords = coords[np.isin(atomnames, ["P", "C1'", "C2"])]
        nucleic_idx = np.where(np.isin(atomnames, ["P", "C1'", "C2"]))[0]
        # Initialize the array to store the indices of the interacting residues
        contacts = []
        # Loop over the protein residues
        for i, protein_coord in enumerate(protein_coords):
            # Compute the distances between the protein residue and all nucleic acid residues
            distances = np.linalg.norm(nucleic_coords - protein_coord, axis=1)
            # if the lowest distance is more than the cutoff, skip
            if np.min(distances) > cutoff:
                continue
            selected_idx = []
            selected_res = []
            for k in np.argsort(distances):
                if distances[k] > cutoff:
                    break
                if (chain_lbl[nucleic_idx[k]],resnums[nucleic_idx[k]]) not in selected_res:
                    selected_idx.append(k)
                    selected_res.append((chain_lbl[nucleic_idx[k]],resnums[nucleic_idx[k]]))
            for k in selected_idx:
                contacts.append([protein_idx[i], nucleic_idx[k]])



            # Find the bead with the minimum distance for each
            # min_distance_idx = np.argmin(distances)
            # Store the indices of the interacting residues
            # contacts.append([protein_idx[i], nucleic_idx[min_distance_idx]])
        return contacts
       
    def construct_hessian(self, pdb_code, build_gamma=False, verbose=False):
        
        if self.coord_dict[pdb_code]:
            coords, chain_lbl, resnames, atomnames, _ = self.coord_dict[pdb_code]

        # print(coords.shape)
        n_atoms = len(coords)
        # print(n_atoms)
        hessian = np.zeros((n_atoms*3, n_atoms*3), float)
        if build_gamma:
            gamma_matrix = np.zeros((n_atoms, n_atoms), float)
        # distance_mat = np.ones((n_atoms*3, n_atoms*3), float)

        #cutoff
        # edENM_cutoff = int(2.9*np.log(n_atoms) - 2.9) # as in Lauras and Domenicos code.
        
        # if edENM_cutoff >= 20:
        #     edENM_cutoff = 20
            
        # if edENM_cutoff <= 8:
        #     edENM_cutoff = 8
        n_protein_atoms = len([c for c in chain_lbl if self.chain_entity_map[c] == "protein"])
        edENM_cutoff = int(np.clip(2.9*np.log(n_protein_atoms) - 2.9, 8, 20))
        edENM_seq_exp = 2 # so we treat neighbours different in sequence.

        n_nucleic_atoms = len([c for c in chain_lbl if self.chain_entity_map[c] == "nucleic"])
        edENM_cutoff_nucleic = 11

        edENM_cutoff_inter = 8 # old approach: 10
        if verbose:
            print("Number of protein atoms: ", n_protein_atoms)
            print("Number of nucleic atoms: ", n_nucleic_atoms)

        edENM_cart_exp = 6 
        edENM_seq_force_constant = 60
        edENM_cart_force_constant = 6
        edENM_M = 3
        
        if n_nucleic_atoms > 0:
            inter_contacts = self._find_protein_nucleic_contacts(pdb_code,edENM_cutoff_inter)
            if verbose:
                print("Number of protein-nucleic spring: ", len(inter_contacts))
            if len(inter_contacts) < 1:
                print("Error! No protein-nucleic contacts found. Check the structure.")
                sys.exit(1)

        for i in range(len(coords)):
            diff = coords[i+1:, :] - coords[i] # diff between x y z of entry n and entry n-1
            squared_diff = diff**2 # squared diff for each x y and z
            # diff = diff between current x y z and  all other x y z later in sequence. 
            
            for j, s_ij in enumerate(squared_diff.sum(1)): # col wise summation   
                # j is the position of the paired bead relative to i (sequential distance is j+1)
                # s_ij is the squared distance between the two beads
                
                # check if protein-protein, protein-nucleic, or nucleic-nucleic and define the spring constant accordingly
                entity_couple = f"{self.chain_entity_map[chain_lbl[i]]}-{self.chain_entity_map[chain_lbl[j+i+1]]}"
                if entity_couple == "protein-protein":
                    # this is the normal edENM
                    if j + 1 <= edENM_M: 
                        gamma = edENM_seq_force_constant/(j+1)**edENM_seq_exp # then we will use this spring constant.
                    elif s_ij <= edENM_cutoff**2:
                        gamma = (edENM_cart_force_constant/np.sqrt(s_ij))**edENM_cart_exp 
                    else:
                        # we are ouside the cutoff, gamma = 0 so we can skip because the hessian is already initialized to 0
                        continue
                elif entity_couple == "nucleic-nucleic":
                    if s_ij <= edENM_cutoff_nucleic**2:
                        couple_res = f"{resnames[i]}-{resnames[j+i+1]}"
                        couple_atom = f"{atomnames[i]}-{atomnames[j+i+1]}"
                        if (couple_atom == "C1'-C2") & (j == 0): # covalent bond
                            gamma = 290 if resnames[i] in ["C","U","T"] else 120 # C,U,T are pirimidines
                        elif (couple_res in ["C-G","G-C"]) & (couple_atom == "C2-C2") & (np.sqrt(s_ij) > 4.1) & (np.sqrt(s_ij) < 4.4): # hbond
                            gamma = 65
                        elif (couple_res == "G-G") & (couple_atom == "C2-C2") & (np.sqrt(s_ij) > 6.7) & (np.sqrt(s_ij) < 7.0): # hbond
                            gamma = 65
                        elif (couple_atom in ["P-P","P-C1'","C1'-P"]) & (chain_lbl[i]==chain_lbl[j+i+1]): # pseudo covalent bond
                            gamma = (20/np.sqrt(s_ij))**2.8
                        else: # vdw
                            gamma = 25/np.sqrt(s_ij)
                    else:
                        # we are ouside the cutoff, gamma = 0 so we can skip because the hessian is already initialized to 0
                        continue
                elif s_ij <= edENM_cutoff_inter**2:
                    protein_idx = i if self.chain_entity_map[chain_lbl[i]] == "protein" else j+i+1
                    nucl_idx = i if self.chain_entity_map[chain_lbl[i]] == "nucleic" else j+i+1
                    if [protein_idx, nucl_idx] not in inter_contacts:
                        continue
                    # protein-nucleic
                    nucl_bead_name = atomnames[nucl_idx]
                    if nucl_bead_name == "C1'":
                        gamma = 30/np.sqrt(s_ij)
                    else:
                        protein_resname = resnames[protein_idx]
                        protein_restype = "CN" if protein_resname in ["ASP","GLU"] else "CP" if protein_resname in ["LYS","ARG","HIS","HIP",'HIE',"HID"] \
                                else "P" if protein_resname in ["SER","THR","ASN","GLN","CYS"] else "H"

                        if ((nucl_bead_name=="C2") & (protein_restype=="CN")) | ((nucl_bead_name=="P") & (protein_restype in ["P","CP"])):
                            gamma =  (25/np.sqrt(s_ij))**2.2 # old approach (10/np.sqrt(s_ij))**1.3
                        else:
                            gamma = 30/np.sqrt(s_ij)
                else:
                    continue
                
                # once the spring constant is defined, we can compute the derivative and update the hessian matrix
                diff_coords = diff[j]
                j = j + i + 1
                derivative = np.outer(diff_coords, diff_coords)*(float(-gamma)/np.sqrt(s_ij)**(2+self.power))  # delta coords * -gamma/dist**0 in our case  in the paper its -gamma * delta diff  / dist**2
                hessian[i*3:i*3+3, j*3:j*3+3] = derivative
                hessian[j*3:j*3+3, i*3:i*3+3] = derivative #symmetry
                hessian[i*3:i*3+3, i*3:i*3+3] = hessian[i*3:i*3+3, i*3:i*3+3] - derivative
                hessian[j*3:j*3+3, j*3:j*3+3] = hessian[j*3:j*3+3, j*3:j*3+3] - derivative
                if build_gamma:
                    gamma_matrix[i, j] = gamma
                    gamma_matrix[j, i] = gamma
        self.hessian_matrices[pdb_code] = hessian
        if build_gamma:
            self.gamma_matrices[pdb_code] = gamma_matrix
        
    

    def add_masses_to_hessian(self, pdb_code):

        if self.hessian_matrices[pdb_code].size > 0: # not empty
            H = self.hessian_matrices[pdb_code]
            masses = self.masses[pdb_code] # masses per residue.
            N_x3 = H.shape[0] # is symmetrical and 3x num_atoms


        MWH = np.zeros_like(H, dtype=np.float64)
        mass_res_expanded = np.repeat(masses, 3) # same mass for x, y , z
        for i in range(N_x3):
            for j in range(N_x3):
                MWH[i, j] = H[i,j] / np.sqrt(mass_res_expanded[i//3] * mass_res_expanded[j//3])
                
            
        self.hessian_weighted_matrices[pdb_code] = MWH
    
    def solve_eigenproblem(self, pdb_code, weighted=False):
        # given that H is symmetrical, we can leverage eigh from numpy.linalg.

        #check symmetry
        if not self.hessian_matrices[pdb_code].size > 0:
            print("Hessian matrix is needed")
            return
    
        # Symmetrize the Hessian matrix
        if weighted:
            if self.hessian_weighted_matrices[pdb_code].size > 0:
                hessian = self.hessian_weighted_matrices[pdb_code]
        else:
            hessian = self.hessian_matrices[pdb_code]
        
        #print(f"{hessian=}")
        #with open("/home/micnag/bioinformatics/domenico_nma/test_weighted_hessian.txt", "w") as fh_out:
        #    np.savetxt(fh_out, hessian, fmt='%7.4f')
        
        #hessian_sym = (hessian + hessian.T) / 2
    
        # Solve the eigenproblem
        eigenvalues, eigenvectors = np.linalg.eigh(hessian)

        num_non_interesting_modes = 6  # Assuming the first 6 modes are not interesting
        eigenvalues = eigenvalues[num_non_interesting_modes:num_non_interesting_modes+self.modes]
        eigenvectors = eigenvectors[:, num_non_interesting_modes:num_non_interesting_modes+self.modes]
    
        # Store the results
        self.eigenvalues[pdb_code] = eigenvalues
        self.eigenvectors[pdb_code] = eigenvectors
        # now lets apply mass weighting as well

    def compute_vectors_for_NMA(self, pdb_code, num_modes=10):
        # we need to compute the vectors for each residue in the structure based on x y z triple pairs in n modes:
        # e.g if dim = 994*3, 10 we need to split first dim into pieces of 3, compute the x y z direction.
        # then we need to grab the corresponding CA residue. Map this vector ontop of the residue. 
        # and we end up with another point that is the endpoint of the vector.

        if self.eigenvectors and self.eigenvalues:
            
            self.new_positions_after_NMA[pdb_code] = {} # empty dict . for each normal mode we store results.

            CA_coords, _, _, _,_ = self.coord_dict[pdb_code] # coords and chain for the struc.


            tmp_dict = {}
            for mode in range(num_modes):
                eigenvecs_mode = self.eigenvectors[pdb_code][:, mode] # gives back the array of shape n_residues*3 , eigenvalues
                
                vecs_reshaped = eigenvecs_mode.reshape((-1, 3)) # shape (N,3) 
                
                new_positions = CA_coords + vecs_reshaped

                tmp_dict[mode] = vecs_reshaped                
                self.new_positions_after_NMA[pdb_code][f"mode_{mode+1}"] = (CA_coords, new_positions)

            self.eigenvecs_reshaped[pdb_code] = tmp_dict
        else:
            # not present eigvecs and eigvals
            print("Eigenvectors or eigenvalues not found for the PDB code provided.")
    
    def write_pymol_script_for_visualization(self, pdb_code, num_modes=10, scale_factor=75):

        # Define parameters for the cone and cylinder (arrow representation)
        cone_radius = 0.5  # Radius of the base of the cone, define this before using it
        cone_length = 0.8  # Length of the cone
        cylinder_radius = 0.2  # Radius of the cylinder, increased for a thicker appearance

        #print(self.new_positions_after_NMA)
        # might be problematic if num modes asked for is < then present in self.new_positions_after_NMA. needsa  check 
        script_location = f"{self.work_dir}/{os.path.basename(pdb_code[:-4])}_normal_mode_{1}_{num_modes}.pml"

        if os.path.exists(script_location):
            os.remove(script_location)

        for mode in range(num_modes):
            
            CA_coords, new_positions = self.new_positions_after_NMA[pdb_code][f"mode_{(mode+1)}"]
            
            with open(script_location, 'a') as script_file:
                script_file.write(f"load {pdb_code}, {os.path.basename(pdb_code)}\n")
                script_file.write(f"hide everything,  {os.path.basename(pdb_code)}\n")
                script_file.write(f"show cartoon,  {os.path.basename(pdb_code)}\n")
        
                for i, (ca, endpoint) in enumerate(zip(CA_coords, new_positions), start=1):
                    scaled_vector = [(e - c) * scale_factor for e, c in zip(endpoint, ca)]
                    scaled_endpoint = [c + v for c, v in zip(ca, scaled_vector)]
        
                    # Begin Python block
                    script_file.write("python\n")
                    script_file.write("from pymol.cgo import CYLINDER, CONE  # Import CGO primitives\n")
                    script_file.write("from pymol import cmd\n")
                    script_file.write("\n")
                    script_file.write("cgo = [\n")
        
                    # Define the cylinder (arrow tail) in yellow
                    script_file.write(f"    CYLINDER, {ca[0]}, {ca[1]}, {ca[2]}, {scaled_endpoint[0]}, {scaled_endpoint[1]}, {scaled_endpoint[2]}, {cylinder_radius}, 1.0, 1.0, 0.0, 1.0, 1.0, 0.0,\n")
        
                    # Define the full cone (arrow head) in yellow with a sharp tip (radius 0.0)
                    cone_tip = [scaled_endpoint[0] + scaled_vector[0] * cone_length,
                                scaled_endpoint[1] + scaled_vector[1] * cone_length,
                                scaled_endpoint[2] + scaled_vector[2] * cone_length]
                    script_file.write(f"    CONE, {scaled_endpoint[0]}, {scaled_endpoint[1]}, {scaled_endpoint[2]}, {cone_tip[0]}, {cone_tip[1]}, {cone_tip[2]}, {cone_radius}, 0.0, 1.0, 1.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0,\n")
        
                    # End CGO definition
                    script_file.write("]\n")
        
                    # Load the CGO object
                    script_file.write(f"cmd.load_cgo(cgo, 'vector_{i}')\n")
                    script_file.write("python end\n\n")
        
                script_file.write(f"zoom {os.path.basename(pdb_code)}\n")
    
        print(f"PyMOL script generated at: {script_location}")


    def animate_NMA(self, template=None, steps=5, both_sides=True, A=150, verbose=False, outpath='.'):
        # go 5 steps in both directions from eq. position.. multiply 1.1 1.2 1.3 1.4 1.5 and minus towards the direction of NM. 
        # overall quality CSO.
        # specific mode.  go 1vs1 inner mode. # abs for range [-1,1] so abs will do scaling between [0,1]

        for mmcif_code in self.coord_dict:
            # coords, chains, resnames, atomnames = self.coord_dict[mmcif_code]  # Extract original coordinates and chain information
            modes_dict = self.eigenvecs_reshaped[mmcif_code]  # Extract eigenvectors for each mode
            if template is None:
                template = mmcif_code
            if template.endswith(".pdb"):
                parser = PDBParser()
            else:
                parser = MMCIFParser()
            original_structure = parser.get_structure("template", template)

            for mode, mode_eigenvecs in modes_dict.items():
                # get only indices for C-alpha atoms and P atoms
                mode_eigenvecs = mode_eigenvecs[(np.array(self.coord_dict[mmcif_code][3]) == "P")|(np.array(self.coord_dict[mmcif_code][3]) == "CA")]

                animated_structure = Structure.Structure(f'animated_mode_{mode}')

                total_steps = steps * 2 if both_sides else steps  # Total steps considering both directions if both_sides is True
                for step in range(-1*steps, steps + 1):
                    if step == 0:
                        continue  # Skip the equilibrium position to avoid duplication

                    ## Calculate the angle for the sine function
                    #angle = math.pi * step / total_steps
                    ## Use the sine of the angle to modulate the factor
                    #factor = math.sin(angle) * A

                    angle = 2 * math.pi * (step / total_steps)   # ONE FULL sine oscillation
                    factor = math.sin(angle) * A

                    new_model = Model.Model(step + steps)  # Adjust model ID to be positive

                    ca_atom_index = 0  # Initialize C-alpha atom index
                    for orig_chain in original_structure.get_chains():
                        new_chain = Chain.Chain(orig_chain.id)

                        for orig_residue in orig_chain.get_residues():
                            for orig_atom in orig_residue:
                                if orig_atom.get_name() in ['CA','P']:  # Filter for C-alpha atoms
                                    if ca_atom_index >= len(mode_eigenvecs):  # Check for index out of bounds
                                        break
                                    eigenvec = mode_eigenvecs[ca_atom_index]
                                    original_coord = np.array(orig_atom.coord)
                                    delta = eigenvec * factor
                                    new_coord = original_coord + delta

                                    new_atom = Atom.Atom(orig_atom.name, new_coord, orig_atom.bfactor, orig_atom.occupancy, orig_atom.altloc, orig_atom.fullname, orig_atom.serial_number, orig_atom.element)
                                    new_residue = Residue.Residue(orig_residue.id, orig_residue.resname, orig_residue.segid)
                                    new_residue.add(new_atom)
                                    new_chain.add(new_residue)

                                    ca_atom_index += 1  # Increment C-alpha atom index

                        new_model.add(new_chain)
                    animated_structure.add(new_model)

                io = PDBIO()
                io.set_structure(animated_structure)
                if (mmcif_code.endswith(".pdb")) | (mmcif_code.endswith(".cif")):
                    output_file_path = f"{outpath}/{mmcif_code[:-4]}_mode_{mode+1}_animated.pdb"
                else:
                    output_file_path = f"{outpath}/{mmcif_code}_mode_{mode+1}_animated.pdb"
                io.save(output_file_path)
                if verbose:
                    print(f"Animated PDB file for mode {mode+1} saved as {output_file_path}")


class classicENM:
    def __init__(self, work_dir=None, script_dir = None,
                 cutoff_inter = 10,
                 gamma = 1,
                 dr = 8,
                 power = 0,
                 modes=10):
       
        self.work_dir = work_dir
        self.script_dir = script_dir
        self.gamma = float(gamma)
        self.dr = float(dr)
        self.power = float(power)
        self.hessian_matrices = {}
        self.coord_dict = {}
        self.eigenvalues = {}
        self.eigenvectors = {}
        self.modes = modes
        self.masses = {}
        self.hessian_weighted_matrices = {}
        self.new_positions_after_NMA = {} # contains the endpoints of the vectors that are used to map the normal modes 
        self.eigenvecs_reshaped = {}
        self.chain_entity_map = {}
        self.gamma_matrices = {}
        self.cutoff_inter = cutoff_inter

    def _change_resname(self,resname):
        if "A" in resname and resname != "GUA" and resname != "URA":
            newname = "A"
        elif "U" in resname and resname != "GUA":
            newname = "U"
            
        elif "G" in resname:
            newname = "G"
        elif "C" in resname:
            newname = "C"
        elif "T" in resname:
            newname = "T"
        else:
            print(f"Error: residue {resname} not recognized")
            sys.exit(1)
        return newname
   
    def extract_ca_coordinates_and_masses(self, pdb, name=None):
        #helper function to parse a pdb file and  

        lst =  [('VAL',99), ('ILE',113), ('LEU',113), ('GLU',129), ('GLN',128),
                        ('ASP',115), ('ASN',114), ('HIS',137), ('TRP',186), ('PHE',147), ('TYR',163), 
                        ('ARG',156), ('LYS',128), ('SER',87), ('THR',101), ('MET',131), ('ALA',71), 
                        ('GLY',57), ('PRO',97), ('CYS',103)]
        
        aa_dict = defaultdict(lambda: 100, lst) # 100 if non canonical
        protein_residues = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
                            "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]
        nucleic_residues = ["DA", "DC", "DG", "DT", "A", "C", "G", "U", "T", "RA", "RC", "RG",
                            "RU", "ADE", "CYT", "GUA", "THY", "URI"]

        if pdb.endswith(".pdb"):
            parser = PDBParser()
        else:
            parser = MMCIFParser()
        
        if name is None:
            name = pdb.split("/")[-1].split(".")[0]
        
        structure = parser.get_structure("protein", pdb) # modify name ?
        
        coords = []
        chain_lbl = []
        masses = []
        resnames = []
        atomnames = []
        
        for model in structure:
            for chain in model:
                # check if chain is protein
                if np.any([residue.get_resname() in protein_residues for residue in chain]):
                    self.chain_entity_map[chain.get_id()] = "protein"
                    for residue in chain:
                        for atom in residue:
                            if atom.get_name() == "CA":
                                atom_coords = atom.get_coord()
                                coords.append(atom_coords) # coord
                                chain_lbl.append(chain.get_id()) # chain
                                res = residue.get_resname()
                                resnames.append(res) # save the resname for protein-nucleic interaction
                                
                                res_weight = np.float64(aa_dict[res]) #* self.mass_conv_u_to_kg # convert to kg from Dalton/Unit

                                masses.append(res_weight) # mass_in_kg
                                atomnames.append("CA")
                elif np.any([residue.get_resname() in nucleic_residues for residue in chain]):
                    self.chain_entity_map[chain.get_id()] = "nucleic"
                    for residue in chain:
                        tmp = {"P": [], "C1'": [], "C2": []} # in this way we are sure to have P-C1-C2 order -> needed for the ENM
                        for atom in residue:
                            if atom.get_name() in tmp.keys():
                                atom_coords = atom.get_coord()
                                tmp[atom.get_name()] = [atom_coords,
                                                        chain.get_id(),
                                                        self._change_resname(residue.get_resname()), # save resname for definition of the spring constant
                                                        np.float64(aa_dict[residue.get_resname()])/3,
                                                        atom.get_name()]
                        for v in tmp.values():
                            if v:
                                coords.append(v[0])
                                chain_lbl.append(v[1])
                                resnames.append(v[2])
                                masses.append(v[3])
                                atomnames.append(v[4])
                else:
                    continue

        self.coord_dict[name] = (np.array(coords), chain_lbl, resnames, atomnames)
        self.masses[name] = masses

    def _find_protein_nucleic_contacts(self, pdb_code,cutoff):
        """
        Function to find protein-nucleic contacts in a structure. In this model, only one spring is
        set between each protein residue and nucleic acid residue, which is between the C-alpha atom
        of the protein residue and the closest bead of the nucleic acid residue.
        The function returns the indices of the couple of interacting residues, shape (n_contacts, 2).
        """
        # Extract the coordinates of the C-alpha atoms of the protein and the P, C1', and C2 atoms of the nucleic acids
        coords, _, _, atomnames = self.coord_dict[pdb_code]
        protein_coords = coords[np.array(atomnames) == "CA"]
        protein_idx = np.where(np.array(atomnames) == "CA")[0]
        nucleic_coords = coords[np.isin(atomnames, ["P", "C1'", "C2"])]
        nucleic_idx = np.where(np.isin(atomnames, ["P", "C1'", "C2"]))[0]
        # Initialize the array to store the indices of the interacting residues
        contacts = []
        # Loop over the protein residues
        for i, protein_coord in enumerate(protein_coords):
            # Compute the distances between the protein residue and all nucleic acid residues
            distances = np.linalg.norm(nucleic_coords - protein_coord, axis=1)
            # if the lowest distance is more than the cutoff, skip
            if np.min(distances) > cutoff:
                continue
            # Find the nucleic acid residue with the minimum distance
            min_distance_idx = np.argmin(distances)
            # Store the indices of the interacting residues
            contacts.append([protein_idx[i], nucleic_idx[min_distance_idx]])
        return contacts
       
    def construct_hessian(self, pdb_code, build_gamma=False, verbose=False,edProt=True):
        
        if self.coord_dict[pdb_code]:
            coords, chain_lbl, resnames, atomnames = self.coord_dict[pdb_code]

        # print(coords.shape)
        n_atoms = len(coords)
        # print(n_atoms)
        hessian = np.zeros((n_atoms*3, n_atoms*3), float)
        if build_gamma:
            gamma_matrix = np.zeros((n_atoms, n_atoms), float)
        # distance_mat = np.ones((n_atoms*3, n_atoms*3), float)
        
        #cutoff
        # edENM_cutoff = int(2.9*np.log(n_atoms) - 2.9) # as in Lauras and Domenicos code.
        
        # if edENM_cutoff >= 20:
        #     edENM_cutoff = 20
            
        # if edENM_cutoff <= 8:
        #     edENM_cutoff = 8
        n_protein_atoms = len([c for c in chain_lbl if self.chain_entity_map[c] == "protein"])
        edENM_cutoff = int(np.clip(2.9*np.log(n_protein_atoms) - 2.9, 8, 20)) if edProt else 12
        edENM_seq_exp = 2 # so we treat neighbours different in sequence.

        n_nucleic_atoms = len([c for c in chain_lbl if self.chain_entity_map[c] == "nucleic"])
        edENM_cutoff_nucleic = 9

        edENM_cutoff_inter = 10
        if verbose:
            print("Number of protein atoms: ", n_protein_atoms)
            print("Number of nucleic atoms: ", n_nucleic_atoms)

        edENM_cart_exp = 6 
        edENM_seq_force_constant = 60
        edENM_cart_force_constant = 6
        edENM_M = 3
        
        if n_nucleic_atoms > 0:
            inter_contacts = self._find_protein_nucleic_contacts(pdb_code,edENM_cutoff_inter)
            if verbose:
                print("Number of protein-nucleic spring: ", len(inter_contacts))
            if len(inter_contacts) < 1:
                print("Error! No protein-nucleic contacts found. Check the structure.")
                sys.exit(1)

        for i in range(len(coords)):
            diff = coords[i+1:, :] - coords[i] # diff between x y z of entry n and entry n-1
            squared_diff = diff**2 # squared diff for each x y and z
            # diff = diff between current x y z and  all other x y z later in sequence. 
            
            for j, s_ij in enumerate(squared_diff.sum(1)): # col wise summation   
                # j is the position of the paired bead relative to i (sequential distance is j+1)
                # s_ij is the squared distance between the two beads
                
                # check if protein-protein, protein-nucleic, or nucleic-nucleic and define the spring constant accordingly
                entity_couple = f"{self.chain_entity_map[chain_lbl[i]]}-{self.chain_entity_map[chain_lbl[j+i+1]]}"
                if entity_couple == "protein-protein":
                    if edProt:
                        # this is the normal edENM
                        if j + 1 <= edENM_M: 
                            gamma = edENM_seq_force_constant/(j+1)**edENM_seq_exp # then we will use this spring constant.
                        elif s_ij <= edENM_cutoff**2:
                            gamma = (edENM_cart_force_constant/np.sqrt(s_ij))**edENM_cart_exp 
                        else:
                            # we are ouside the cutoff, gamma = 0 so we can skip because the hessian is already initialized to 0
                            continue
                    else:
                        if s_ij <= edENM_cutoff**2:
                            gamma = 1
                        else:
                            continue
                elif (entity_couple == "nucleic-nucleic") & (s_ij <= edENM_cutoff_nucleic**2):
                    gamma = 1
                elif s_ij <= edENM_cutoff_inter**2:
                    gamma = 1
                else:
                    continue
                
                # once the spring constant is defined, we can compute the derivative and update the hessian matrix
                diff_coords = diff[j]
                j = j + i + 1
                derivative = np.outer(diff_coords, diff_coords)*(float(-gamma)/np.sqrt(s_ij)**(2+self.power))  # delta coords * -gamma/dist**0 in our case  in the paper its -gamma * delta diff  / dist**2
                hessian[i*3:i*3+3, j*3:j*3+3] = derivative
                hessian[j*3:j*3+3, i*3:i*3+3] = derivative #symmetry
                hessian[i*3:i*3+3, i*3:i*3+3] = hessian[i*3:i*3+3, i*3:i*3+3] - derivative
                hessian[j*3:j*3+3, j*3:j*3+3] = hessian[j*3:j*3+3, j*3:j*3+3] - derivative
                if build_gamma:
                    gamma_matrix[i, j] = gamma
                    gamma_matrix[j, i] = gamma
        self.hessian_matrices[pdb_code] = hessian
        if build_gamma:
            self.gamma_matrices[pdb_code] = gamma_matrix
        
    

    def add_masses_to_hessian(self, pdb_code):

        if self.hessian_matrices[pdb_code].size > 0: # not empty
            H = self.hessian_matrices[pdb_code]
            masses = self.masses[pdb_code] # masses per residue.
            N_x3 = H.shape[0] # is symmetrical and 3x num_atoms


        MWH = np.zeros_like(H, dtype=np.float64)
        mass_res_expanded = np.repeat(masses, 3) # same mass for x, y , z
        for i in range(N_x3):
            for j in range(N_x3):
                MWH[i, j] = H[i,j] / np.sqrt(mass_res_expanded[i//3] * mass_res_expanded[j//3])
                
            
        self.hessian_weighted_matrices[pdb_code] = MWH
    
    def solve_eigenproblem(self, pdb_code, weighted=False):
        # given that H is symmetrical, we can leverage eigh from numpy.linalg.

        #check symmetry
        if not self.hessian_matrices[pdb_code].size > 0:
            print("Hessian matrix is needed")
            return
    
        # Symmetrize the Hessian matrix
        if weighted:
            if self.hessian_weighted_matrices[pdb_code].size > 0:
                hessian = self.hessian_weighted_matrices[pdb_code]
        else:
            hessian = self.hessian_matrices[pdb_code]
        
        #print(f"{hessian=}")
        #with open("/home/micnag/bioinformatics/domenico_nma/test_weighted_hessian.txt", "w") as fh_out:
        #    np.savetxt(fh_out, hessian, fmt='%7.4f')
        
        #hessian_sym = (hessian + hessian.T) / 2
    
        # Solve the eigenproblem
        eigenvalues, eigenvectors = np.linalg.eigh(hessian)

        num_non_interesting_modes = 6  # Assuming the first 6 modes are not interesting
        eigenvalues = eigenvalues[num_non_interesting_modes:num_non_interesting_modes+self.modes]
        eigenvectors = eigenvectors[:, num_non_interesting_modes:num_non_interesting_modes+self.modes]
    
        # Store the results
        self.eigenvalues[pdb_code] = eigenvalues
        self.eigenvectors[pdb_code] = eigenvectors
        # now lets apply mass weighting as well

    def compute_vectors_for_NMA(self, pdb_code, num_modes=10):
        # we need to compute the vectors for each residue in the structure based on x y z triple pairs in n modes:
        # e.g if dim = 994*3, 10 we need to split first dim into pieces of 3, compute the x y z direction.
        # then we need to grab the corresponding CA residue. Map this vector ontop of the residue. 
        # and we end up with another point that is the endpoint of the vector.

        if self.eigenvectors and self.eigenvalues:
            
            self.new_positions_after_NMA[pdb_code] = {} # empty dict . for each normal mode we store results.

            CA_coords, _, _, _ = self.coord_dict[pdb_code] # coords and chain for the struc.


            tmp_dict = {}
            for mode in range(num_modes):
                eigenvecs_mode = self.eigenvectors[pdb_code][:, mode] # gives back the array of shape n_residues*3 , eigenvalues
                
                vecs_reshaped = eigenvecs_mode.reshape((-1, 3)) # shape (N,3) 
                
                new_positions = CA_coords + vecs_reshaped

                tmp_dict[mode] = vecs_reshaped                
                self.new_positions_after_NMA[pdb_code][f"mode_{mode+1}"] = (CA_coords, new_positions)

            self.eigenvecs_reshaped[pdb_code] = tmp_dict
        else:
            # not present eigvecs and eigvals
            print("Eigenvectors or eigenvalues not found for the PDB code provided.")
    
    def write_pymol_script_for_visualization(self, pdb_code, num_modes=10, scale_factor=75):

        # Define parameters for the cone and cylinder (arrow representation)
        cone_radius = 0.5  # Radius of the base of the cone, define this before using it
        cone_length = 0.8  # Length of the cone
        cylinder_radius = 0.2  # Radius of the cylinder, increased for a thicker appearance

        #print(self.new_positions_after_NMA)
        # might be problematic if num modes asked for is < then present in self.new_positions_after_NMA. needsa  check 
        script_location = f"{self.work_dir}/{os.path.basename(pdb_code[:-4])}_normal_mode_{1}_{num_modes}.pml"

        if os.path.exists(script_location):
            os.remove(script_location)

        for mode in range(num_modes):
            
            CA_coords, new_positions = self.new_positions_after_NMA[pdb_code][f"mode_{(mode+1)}"]
            
            with open(script_location, 'a') as script_file:
                script_file.write(f"load {pdb_code}, {os.path.basename(pdb_code)}\n")
                script_file.write(f"hide everything,  {os.path.basename(pdb_code)}\n")
                script_file.write(f"show cartoon,  {os.path.basename(pdb_code)}\n")
        
                for i, (ca, endpoint) in enumerate(zip(CA_coords, new_positions), start=1):
                    scaled_vector = [(e - c) * scale_factor for e, c in zip(endpoint, ca)]
                    scaled_endpoint = [c + v for c, v in zip(ca, scaled_vector)]
        
                    # Begin Python block
                    script_file.write("python\n")
                    script_file.write("from pymol.cgo import CYLINDER, CONE  # Import CGO primitives\n")
                    script_file.write("from pymol import cmd\n")
                    script_file.write("\n")
                    script_file.write("cgo = [\n")
        
                    # Define the cylinder (arrow tail) in yellow
                    script_file.write(f"    CYLINDER, {ca[0]}, {ca[1]}, {ca[2]}, {scaled_endpoint[0]}, {scaled_endpoint[1]}, {scaled_endpoint[2]}, {cylinder_radius}, 1.0, 1.0, 0.0, 1.0, 1.0, 0.0,\n")
        
                    # Define the full cone (arrow head) in yellow with a sharp tip (radius 0.0)
                    cone_tip = [scaled_endpoint[0] + scaled_vector[0] * cone_length,
                                scaled_endpoint[1] + scaled_vector[1] * cone_length,
                                scaled_endpoint[2] + scaled_vector[2] * cone_length]
                    script_file.write(f"    CONE, {scaled_endpoint[0]}, {scaled_endpoint[1]}, {scaled_endpoint[2]}, {cone_tip[0]}, {cone_tip[1]}, {cone_tip[2]}, {cone_radius}, 0.0, 1.0, 1.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0,\n")
        
                    # End CGO definition
                    script_file.write("]\n")
        
                    # Load the CGO object
                    script_file.write(f"cmd.load_cgo(cgo, 'vector_{i}')\n")
                    script_file.write("python end\n\n")
        
                script_file.write(f"zoom {os.path.basename(pdb_code)}\n")
    
        print(f"PyMOL script generated at: {script_location}")


    def animate_NMA(self, template=None, steps=5, both_sides=True, A=150, verbose=False, outpath='.'):
        # go 5 steps in both directions from eq. position.. multiply 1.1 1.2 1.3 1.4 1.5 and minus towards the direction of NM. 
        # overall quality CSO.
        # specific mode.  go 1vs1 inner mode. # abs for range [-1,1] so abs will do scaling between [0,1]

        for mmcif_code in self.coord_dict:
            # coords, chains, resnames, atomnames = self.coord_dict[mmcif_code]  # Extract original coordinates and chain information
            modes_dict = self.eigenvecs_reshaped[mmcif_code]  # Extract eigenvectors for each mode
            if template is None:
                template = mmcif_code
            if template.endswith(".pdb"):
                parser = PDBParser()
            else:
                parser = MMCIFParser()
            original_structure = parser.get_structure("template", template)

            for mode, mode_eigenvecs in modes_dict.items():
                # get only indices for C-alpha atoms and P atoms
                mode_eigenvecs = mode_eigenvecs[(np.array(self.coord_dict[mmcif_code][3]) == "P")|(np.array(self.coord_dict[mmcif_code][3]) == "CA")]

                animated_structure = Structure.Structure(f'animated_mode_{mode}')

                total_steps = steps * 2 if both_sides else steps  # Total steps considering both directions if both_sides is True
                for step in range(-1*steps, steps + 1):
                    if step == 0:
                        continue  # Skip the equilibrium position to avoid duplication

                    ## Calculate the angle for the sine function
                    #angle = math.pi * step / total_steps
                    ## Use the sine of the angle to modulate the factor
                    #factor = math.sin(angle) * A

                    angle = 2 * math.pi * (step / total_steps)
                    factor = math.sin(angle) * A

                    new_model = Model.Model(step + steps)  # Adjust model ID to be positive

                    ca_atom_index = 0  # Initialize C-alpha atom index
                    for orig_chain in original_structure.get_chains():
                        new_chain = Chain.Chain(orig_chain.id)

                        for orig_residue in orig_chain.get_residues():
                            for orig_atom in orig_residue:
                                if orig_atom.get_name() in ['CA','P']:  # Filter for C-alpha atoms
                                    if ca_atom_index >= len(mode_eigenvecs):  # Check for index out of bounds
                                        break
                                    eigenvec = mode_eigenvecs[ca_atom_index]
                                    original_coord = np.array(orig_atom.coord)
                                    delta = eigenvec * factor
                                    new_coord = original_coord + delta

                                    new_atom = Atom.Atom(orig_atom.name, new_coord, orig_atom.bfactor, orig_atom.occupancy, orig_atom.altloc, orig_atom.fullname, orig_atom.serial_number, orig_atom.element)
                                    new_residue = Residue.Residue(orig_residue.id, orig_residue.resname, orig_residue.segid)
                                    new_residue.add(new_atom)
                                    new_chain.add(new_residue)

                                    ca_atom_index += 1  # Increment C-alpha atom index

                        new_model.add(new_chain)
                    animated_structure.add(new_model)

                io = PDBIO()
                io.set_structure(animated_structure)
                if (mmcif_code.endswith(".pdb")) | (mmcif_code.endswith(".cif")):
                    output_file_path = f"{outpath}/{mmcif_code[:-4]}_mode_{mode+1}_animated.pdb"
                else:
                    output_file_path = f"{outpath}/{mmcif_code}_mode_{mode+1}_animated.pdb"
                io.save(output_file_path)
                if verbose:
                    print(f"Animated PDB file for mode {mode+1} saved as {output_file_path}")