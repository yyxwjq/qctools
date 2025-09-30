'''
Author: glab-cabage 2227541807@qq.com
Date: 2024-04-26 02:52:52
LastEditors: glab-cabage 2227541807@qq.com
LastEditTime: 2024-07-24 16:01:31
'''

import os
import logging
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
try:
    from pynep.calculate import NEP
except ImportError:
    print("If you want to use nep calculator,please install pynep first!")
import matplotlib.ticker as ticker
from ase.io import read, write, Trajectory
from matplotlib.ticker import MultipleLocator
from qctools import qctools_logging

# Initialize logging
qctools_logging()

'''
apps: n2p2 or nep                                                type: str
data: images data or training or testing results from software.  type: class<Atoms> or file np.loadtxt can read.
(depends on resource parameters)
resource: software or images                                     type: str
'''

def data_normalization(apps, 
                       imgdata, 
                       resource, 
                       data=None, 
                       pot=None
):

# NOTE: Currently only supports n2p2 and nep.
    if resource == 'software':
        imgs = imgdata
        if data is not None:
            # Support both string and dict format for data parameter
            if isinstance(data, dict):
                # Dictionary format: {'energy': 'energy_file.txt', 'force': 'force_file.txt'}
                energy_file = data.get('energy')
                force_file = data.get('force')
                
                # Process energy data if provided
                if energy_file:
                    try:
                        e_dat = np.loadtxt(energy_file)
                        logging.info(f'Successfully loaded energy data: {energy_file}, shape: {e_dat.shape}')
                        
                        if apps == 'nep':
                            if e_dat.shape[1] == 2 and e_dat.shape[0] == len(imgs):
                                img_num_matrix = np.arange(0, len(imgs)).reshape(-1, 1)
                                normalized_array = np.concatenate((img_num_matrix, e_dat[:, [1, 0]]), axis=1) # [1, 0] DFT NEP
                                np.savetxt('energy.data', normalized_array, fmt='%4s')
                                logging.info(f'Energy data normalized and saved to energy.data ({len(imgs)} structures)')
                            else:
                                logging.error(f'Energy data format error: expected shape ({len(imgs)}, 2), got {e_dat.shape}')
                        elif apps == 'n2p2':
                            if e_dat.shape[1] == 3:
                                np.savetxt('energy.data', e_dat, fmt='%4s')
                                logging.info('n2p2 energy data saved to energy.data')
                            else:
                                logging.error(f'n2p2 energy data format error: expected 3 columns, got {e_dat.shape[1]}')
                    except Exception as e:
                        logging.error(f'Failed to load energy data from {energy_file}: {e}')
                
                # Process force data if provided
                if force_file:
                    try:
                        f_dat = np.loadtxt(force_file)
                        logging.info(f'Successfully loaded force data: {force_file}, shape: {f_dat.shape}')
                        
                        if apps == 'nep':
                            atom_num = sum(map(lambda x : len(x), imgs))
                            if f_dat.shape[1] == 6 and f_dat.shape[0] == atom_num:
                                img_atom_matrix = []
                                for i in range(len(imgs)):
                                    for j in range(len(imgs[i])):
                                        img_atom_matrix.append([i, j])
                                normalized_array = np.concatenate((img_atom_matrix, f_dat[:, [3, 4, 5, 0, 1, 2]]), axis=1)
                                np.savetxt('force.data', normalized_array, fmt='%4s')
                                logging.info(f'Force data normalized and saved to force.data ({atom_num} atoms)')
                            else:
                                logging.error(f'Force data format error: expected shape ({atom_num}, 6), got {f_dat.shape}')
                        elif apps == 'n2p2':
                            if f_dat.shape[1] == 4:
                                atom_num = int(f_dat.shape[0] / 3)
                                normalized_array = []
                                for i in range(atom_num):
                                    line = [f_dat[i*3+2][0], f_dat[i*3+2][1], f_dat[i*3][2], f_dat[i*3+1][2], f_dat[i*3+2][2], f_dat[i*3][3], f_dat[i*3+1][3], f_dat[i*3+2][3]]
                                    normalized_array.append(line)
                                np.savetxt('force.data', normalized_array, fmt='%4s')
                                logging.info(f'n2p2 force data normalized and saved to force.data ({atom_num} atoms)')
                            else:
                                logging.error(f'n2p2 force data format error: expected 4 columns, got {f_dat.shape[1]}')
                    except Exception as e:
                        logging.error(f'Failed to load force data from {force_file}: {e}')
                        
            else:
                # Legacy string format - try to determine if it's energy or force based on data shape
                try:
                    e_f_dat = np.loadtxt(data)
                    logging.info(f'Successfully loaded training data: {data}, shape: {e_f_dat.shape}')
                    
                    if apps == 'nep':
                        if e_f_dat.shape[1] == 2:
                            if e_f_dat.shape[0] == len(imgs):
                                img_num_matrix = np.arange(0, len(imgs)).reshape(-1, 1)
                                normalized_array = np.concatenate((img_num_matrix, e_f_dat[:, [1, 0]]), axis=1) # [1, 0] DFT NEP
                                np.savetxt('energy.data', normalized_array, fmt='%4s')
                                logging.info(f'Energy data normalized and saved to energy.data ({len(imgs)} structures)')
                            else:
                                logging.error(f'Data count mismatch: {e_f_dat.shape[0]} vs {len(imgs)}')
                                return
                                
                        elif e_f_dat.shape[1] == 6:
                            atom_num = sum(map(lambda x : len(x), imgs))
                            if e_f_dat.shape[0] == atom_num:
                                img_atom_matrix = []
                                for i in range(len(imgs)):
                                    for j in range(len(imgs[i])):
                                        img_atom_matrix.append([i, j])
                                normalized_array = np.concatenate((img_atom_matrix, e_f_dat[:, [3, 4, 5, 0, 1, 2]]), axis=1)
                                np.savetxt('force.data', normalized_array, fmt='%4s')
                                logging.info(f'Force data normalized and saved to force.data ({atom_num} atoms)')
                            else:
                                logging.error(f'Atom count mismatch: {e_f_dat.shape[0]} vs {atom_num}')
                                return
                                
                    elif apps == 'n2p2':
                        if e_f_dat.shape[1] == 3:
                            np.savetxt('energy.data', e_f_dat, fmt='%4s')
                            logging.info('n2p2 energy data saved to energy.data')
                        elif e_f_dat.shape[1] == 4:
                            atom_num = int(e_f_dat.shape[0] / 3)
                            normalized_array = []
                            for i in range(atom_num):
                                line = [e_f_dat[i*3+2][0], e_f_dat[i*3+2][1], e_f_dat[i*3][2], e_f_dat[i*3+1][2], e_f_dat[i*3+2][2], e_f_dat[i*3][3], e_f_dat[i*3+1][3], e_f_dat[i*3+2][3]]
                                normalized_array.append(line)
                            
                            np.savetxt('force.data', normalized_array, fmt='%4s')
                            logging.info(f'n2p2 force data normalized and saved to force.data ({atom_num} atoms)')
                    else:
                        logging.error(f'Unsupported ML software: {apps}')
                        
                except Exception as e:
                    logging.error(f'Failed to load training data: {e}')
                    return
        else:
            logging.error(f'Need to load training results file for {apps}')
            return
                
# TODO: When the input are structures. How to generate the normalized data.                  
    elif resource == 'images':
        imgs = imgdata
        logging.info(f'Processing {len(imgs)} structures using ML potential: {pot}')
        if apps == 'nep':
            calc = NEP(pot)
        e_dft, e_mlp, f_dft, f_mlp = [], [], [], []
        img_num_matrix = np.arange(0, len(imgs)).reshape(-1, 1)
        img_atom_matrix = []
        for i, img in enumerate(imgs):
            img_copy = img.copy()
            e_dft.append(img.get_potential_energy() / len(img))
            f_dft.append(img.get_forces())
            img_copy.calc = calc
            e_mlp.append(img_copy.get_potential_energy() / len(img))
            f_mlp.append(img_copy.get_forces())
            for j in range(len(img)):
                img_atom_matrix.append([i, j])
        e_dft = np.array(e_dft).reshape(-1, 1)
        e_mlp = np.array(e_mlp).reshape(-1, 1)  # Fixed bug: was e_dft
        f_dft = np.concatenate(f_dft, axis=0)
        f_mlp = np.concatenate(f_mlp, axis=0)
        img_atom_matrix = np.array(img_atom_matrix)
        # print(np.shape(e_dft), np.shape(e_mlp), np.shape(img_num_matrix), np.shape(f_dft), np.shape(f_mlp), np.shape(img_atom_matrix))
        ene_normalized_array = np.concatenate((img_num_matrix, e_dft, e_mlp), axis=1)
        force_normalized_array = np.concatenate((img_atom_matrix, f_dft, f_mlp), axis=1)
        logging.info(f'Data shapes - Energy: {ene_normalized_array.shape}, Force: {force_normalized_array.shape}')
        np.savetxt('energy.data', ene_normalized_array, fmt='%7s')
        np.savetxt('force.data', force_normalized_array, fmt='%7s')
        logging.info('Data saved to energy.data and force.data')

def rmse(T, P):
    T, P = np.array(T), np.array(P)
    error = T - P
    RMSE = np.sqrt(np.mean(np.square(error)))
    return error, RMSE

def err_structure_finding(error_bar, 
                          images, 
                          fontsize, 
                          replace_atom, 
                          Cutimg=False, 
                          comment=False,
                          show_marginals=True
):
    
    # Process energy data
    energy_data = None
    try:
        energy_data = np.loadtxt('energy.data')
        logging.info(f'Loaded energy data: {energy_data.shape}')
        
        error, RMSE = rmse(energy_data[:, 1], energy_data[:, 2])
        logging.info(f'Energy RMSE: {RMSE:.6f} eV/atom')
        err_indices = np.array(np.where(np.abs(error) > error_bar * RMSE)[0])
        logging.info(f'Found {len(err_indices)} structures with large energy errors (>{error_bar} × RMSE)')
        
        if len(err_indices) > 0:
            np.savetxt('energy_error.txt', energy_data[err_indices], fmt='%4s')
            err_img = [images[index] for index in err_indices]
            try:
                write('Err-energy.xyz', err_img, format='extxyz')
            except Exception as e:
                logging.error(f"Error writing Err-energy.xyz: {e}")
        else:
            logging.info('No structures with large energy errors found.')
        
        if Cutimg:
            leave_img_id = [x for x in range(len(images)) if x not in err_indices]
            leave_img = [images[index] for index in leave_img_id]
            write('leave-E-img.xyz', leave_img, format='extxyz')
            
    except Exception as e:
        logging.error(f"Error loading energy.data: {e}")
    
    # Process force data
    force_data = None
    try:
        force_data = np.loadtxt('force.data')
        logging.info(f'Loaded force data: {force_data.shape}')
        
        error, RMSE = rmse(force_data[:, 2:5], force_data[:, 5:8])
        logging.info(f'Force RMSE: {RMSE:.6f} eV/Å')
        frr_indices = np.argwhere(np.abs(error) > error_bar * RMSE)
        logging.info(f'Found {len(frr_indices)} force components with large errors (>{error_bar} × RMSE)')
        
        if len(frr_indices) > 0:
            frr_output = [[force_data[xx[0], 0], force_data[xx[0], 1], force_data[xx[0], xx[1]+2], force_data[xx[0], xx[1]+5]] for xx in frr_indices]
            np.savetxt('force_error.txt', frr_output, fmt='%4s')
            
            bank = np.array(frr_output)[:, 0:2]
            bank1 = np.unique(bank, axis=0).astype(float).astype(int)
            dic = {}
            for i in bank1:
                if i[0] in dic:
                    dic[i[0]].append(i[1])
                else:
                    dic[i[0]] = [i[1]]
            
            img_info = sorted(dic.items())
            
            # Remove existing error structure files
            for filename in ['Err-force-ini.xyz', 'Err-force-replaced.xyz']:
                if os.path.exists(filename):
                    os.remove(filename)

            # Prepare original and replaced structure lists
            original_structures = []
            replaced_structures = []
            
            logging.info(f'Processing {len(img_info)} error structures...')
            
            for k in img_info:
                try:
                    # Get original structure
                    original_structure = images[k[0]]
                    original_structures.append(original_structure)
                    
                    # Create replaced structure
                    from ase import Atoms
                    positions = original_structure.get_positions()
                    symbols = original_structure.get_chemical_symbols()
                    cell = original_structure.get_cell()
                    pbc = original_structure.get_pbc()
                    
                    # Create replaced symbols
                    replaced_symbols = symbols.copy()
                    for j in k[1]:
                        try:
                            if j >= len(replaced_symbols):
                                logging.warning(f'Atom index {j} out of range for structure {k[0]} (has {len(replaced_symbols)} atoms)')
                                continue
                                
                            if isinstance(replace_atom, str):
                                replaced_symbols[j] = replace_atom
                            elif isinstance(replace_atom, dict):
                                original_symbol = symbols[j]
                                if original_symbol in replace_atom:
                                    replaced_symbols[j] = replace_atom[original_symbol]
                                else:
                                    logging.warning(f'No replacement defined for atom type {original_symbol}')
                            else:
                                logging.error(f'Invalid replace_atom parameter type: {type(replace_atom)}')
                                break
                        except (IndexError, KeyError) as e:
                            logging.error(f'Error replacing atom (structure {k[0]}, position {j}): {e}')
                            continue
                    
                    # Create replaced structure
                    replaced_structure = Atoms(
                        symbols=replaced_symbols,
                        positions=positions.copy(),
                        cell=cell.copy(),
                        pbc=pbc.copy()
                    )
                    
                    # Copy additional properties carefully
                    if hasattr(original_structure, 'info') and original_structure.info:
                        replaced_structure.info.update(original_structure.info)
                    
                    replaced_structures.append(replaced_structure)
                    
                except Exception as e:
                    logging.error(f'Error processing structure {k[0]}: {e}')
                    continue
            
            # Write structures using safe function
            write('Err-force-ini.xyz', original_structures, format='extxyz')
            write('Err-force-replaced.xyz', replaced_structures, format='extxyz')
            
            err_img = [info[0] for info in img_info]
            if Cutimg:
                leave_img_id = [x for x in range(len(images)) if x not in err_img]
                leave_img = [images[index] for index in leave_img_id]
                write('leave-F-img.xyz', leave_img, format='extxyz')
        else:
            logging.info('No structures with large force errors found.')
            
    except Exception as e:
        logging.error(f"Error loading force.data: {e}")
    
    # Generate plots for both energy and force if data is available
    if energy_data is not None:
        plot_scatter_with_marginals(x=energy_data[:, 1].reshape(-1, 1), 
                                  y=energy_data[:, 2].reshape(-1, 1), 
                                  fontsize=fontsize, 
                                  plot_type='energy',
                                  Err_comment=comment,
                                  show_marginals=show_marginals,
                                  colors=['#ef5674', '#eebc59','#5891d5'])
    
    if force_data is not None:
        plot_scatter_with_marginals(x=force_data[:, 2:5], 
                                  y=force_data[:, 5:8], 
                                  fontsize=fontsize, 
                                  plot_type='force',
                                  Err_comment=comment,
                                  show_marginals=show_marginals,
                                  colors=['#ef5674', '#eebc59','#5891d5'])
    
    
def plot_scatter_with_marginals(x, 
                               y, 
                               fontsize, 
                               plot_type='energy',
                               Err_comment=True, 
                               show_marginals=True,
                               colors=['blue','green','orange'], 
):
    """Plot scatter with optional marginal density distributions"""
    
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec
    import matplotlib.ticker as ticker
    from scipy import stats
    
    # Base settings
    matplotlib.use('Agg')
    font3 = {'family':'Times','size':fontsize,'color':'k'}
    matplotlib.rcParams['xtick.direction'] = 'in'
    matplotlib.rcParams['ytick.direction'] = 'in'
    matplotlib.rcParams['ytick.right'] = True
    matplotlib.rcParams['xtick.top'] = True
    plt.rcParams['font.size'] = fontsize
    plt.rcParams['font.family'] = 'Times'
    plt.rcParams['font.weight'] = "heavy"
    
    # Create figure with different layouts based on show_marginals
    if show_marginals:
        fig = plt.figure(figsize=(8, 8))
        gs = GridSpec(4, 4, figure=fig, hspace=0.0, wspace=0.0)
        
        # Main scatter plot
        ax_main = fig.add_subplot(gs[1:, :-1])
        # Top histogram (x-axis marginal)
        ax_top = fig.add_subplot(gs[0, :-1], sharex=ax_main)
        # Right histogram (y-axis marginal)
        ax_right = fig.add_subplot(gs[1:, -1], sharey=ax_main)
    else:
        fig, ax_main = plt.subplots(figsize=(8, 8))
        ax_top = None
        ax_right = None
    
    # Calculate RMSE and limits
    error, RMSE = rmse(x, y)
    lowLimit = np.min([np.min(x), np.min(y)])
    highLimit = np.max([np.max(x), np.max(y)])
    
    # Main scatter plot
    if x.shape[1] == y.shape[1] == 1:
        ax_main.scatter(x, y, c=colors[2], alpha=0.8, s=15)
        ax_main.set_ylabel("DFT energy (eV/atom)", fontdict=font3)
        ax_main.set_xlabel("NEP energy (eV/atom)", fontdict=font3)
        ax_main.text(highLimit, lowLimit,
                    "Energy RMSE: {:.3f} eV/atom".format(RMSE),
                    ha='right', va='bottom', fontsize=fontsize-2)
        
        # Marginal density plots for energy
        if show_marginals:
            # Top marginal (x-axis)
            x_flat = x.flatten()
            x_range = np.linspace(x_flat.min(), x_flat.max(), 100)
            try:
                kde_x = stats.gaussian_kde(x_flat)
                density_x = kde_x(x_range)
                ax_top.fill_between(x_range, density_x, alpha=0.4, color=colors[0])
                ax_top.plot(x_range, density_x, color=colors[0], linewidth=2, alpha=0.9)
            except Exception as e:
                logging.warning(f'Failed to create KDE for x-axis: {e}, using histogram instead')
                ax_top.hist(x_flat, bins=30, alpha=0.7, color=colors[0], density=True)
            
            # Right marginal (y-axis)
            y_flat = y.flatten()
            y_range = np.linspace(y_flat.min(), y_flat.max(), 100)
            try:
                kde_y = stats.gaussian_kde(y_flat)
                density_y = kde_y(y_range)
                ax_right.fill_betweenx(y_range, density_y, alpha=0.4, color=colors[1])
                ax_right.plot(density_y, y_range, color=colors[1], linewidth=2, alpha=0.9)
            except Exception as e:
                logging.warning(f'Failed to create KDE for y-axis: {e}, using histogram instead')
                ax_right.hist(y_flat, bins=30, alpha=0.7, color=colors[1], density=True, orientation='horizontal')
        
    elif x.shape[1] == y.shape[1] == 3:
        # Flatten arrays for force components
        x_flat = x.flatten()
        y_flat = y.flatten()
        lowLimit = np.min([np.min(x_flat), np.min(y_flat)])
        highLimit = np.max([np.max(x_flat), np.max(y_flat)])
        
        # Create color map for force components
        color_map = [colors[i % 3] for i in range(3) for _ in range(x.shape[0])]
        
        # Create scatter plot with labels for legend
        for i in range(3):
            start_idx = i * x.shape[0]
            end_idx = (i + 1) * x.shape[0]
            ax_main.scatter(x_flat[start_idx:end_idx], y_flat[start_idx:end_idx], 
                           c=colors[i], alpha=0.8, s=12, 
                           label=f'F_{["x", "y", "z"][i]}')
        
        # Add legend for force components
        ax_main.legend(loc='upper left', fontsize=fontsize-4, framealpha=0.8)
        ax_main.set_ylabel("DFT force (eV/Å)", fontdict=font3)
        ax_main.set_xlabel("NEP force (eV/Å)", fontdict=font3)
        ax_main.text(highLimit, lowLimit,
                    "Force RMSE: {:.3f} eV/Å".format(RMSE),
                    ha='right', va='bottom', fontsize=fontsize-2)
        
        # Marginal density plots for force (show all three components)
        if show_marginals:
            # Top marginal (x-axis) - show all three force components
            x_range = np.linspace(x_flat.min(), x_flat.max(), 100)
            try:
                # Plot density for each force component separately
                for i in range(3):
                    component_data = x[:, i]  # Get x, y, z components separately
                    kde_x = stats.gaussian_kde(component_data)
                    density_x = kde_x(x_range)
                    ax_top.fill_between(x_range, density_x, alpha=0.3, color=colors[i], 
                                       label=f'F_{["x", "y", "z"][i]}')
                    ax_top.plot(x_range, density_x, color=colors[i], linewidth=2, alpha=0.9)
            except Exception as e:
                logging.warning(f'Failed to create KDE for x-axis: {e}, using histogram instead')
                # Fallback to histogram if KDE fails
                for i in range(3):
                    component_data = x[:, i]
                    ax_top.hist(component_data, bins=30, alpha=0.4, color=colors[i], 
                              density=True, histtype='step', linewidth=2)
            
            # Right marginal (y-axis) - show all three force components
            y_range = np.linspace(y_flat.min(), y_flat.max(), 100)
            try:
                # Plot density for each force component separately
                for i in range(3):
                    component_data = y[:, i]  # Get x, y, z components separately
                    kde_y = stats.gaussian_kde(component_data)
                    density_y = kde_y(y_range)
                    ax_right.fill_betweenx(y_range, density_y, alpha=0.3, color=colors[i])
                    ax_right.plot(density_y, y_range, color=colors[i], linewidth=2, alpha=0.9)
            except Exception as e:
                logging.warning(f'Failed to create KDE for y-axis: {e}, using histogram instead')
                # Fallback to histogram if KDE fails
                for i in range(3):
                    component_data = y[:, i]
                    ax_right.hist(component_data, bins=30, alpha=0.4, color=colors[i], 
                                density=True, orientation='horizontal', histtype='step', linewidth=2)
        
    else:
        logging.warning('Using default X-Y label settings')
        ax_main.scatter(x, y, c=colors[2], alpha=0.7, s=20)
        if show_marginals:
            x_flat = x.flatten()
            y_flat = y.flatten()
            # Simple fallback for unknown data format
            ax_top.hist(x_flat, bins=30, alpha=0.7, color=colors[0], density=True)
            ax_right.hist(y_flat, bins=30, alpha=0.7, color=colors[1], density=True, orientation='horizontal')
    
    # Draw diagonal line
    ax_main.plot([lowLimit, highLimit], [lowLimit, highLimit], '--', c='black', linewidth=1)
    
    # Format main plot
    xmajorLocator = ticker.MaxNLocator(5)
    ymajorLocator = ticker.MaxNLocator(5)
    ax_main.xaxis.set_major_locator(xmajorLocator)
    ax_main.yaxis.set_major_locator(ymajorLocator)
    
    xmajorFormatter = ticker.FormatStrFormatter('%.2f')
    ymajorFormatter = ticker.FormatStrFormatter('%.2f')
    ax_main.xaxis.set_major_formatter(xmajorFormatter)
    ax_main.yaxis.set_major_formatter(ymajorFormatter)
    
    for spine in ax_main.spines.values():
        spine.set_linewidth(2)
    
    ax_main.tick_params(labelsize=fontsize-2)
    for label in ax_main.get_xticklabels() + ax_main.get_yticklabels():
        label.set_fontweight('heavy')
    
    # Format marginal plots if they exist
    if show_marginals and ax_top is not None and ax_right is not None:
        # Remove tick labels for marginal plots
        ax_top.tick_params(labelbottom=False, labeltop=False, labelright=False, labelleft=False)
        ax_right.tick_params(labelbottom=False, labeltop=False, labelright=False, labelleft=False)
        
        # Remove all spines except the ones connecting to main plot
        ax_top.spines['top'].set_visible(False)
        ax_top.spines['right'].set_visible(False)
        ax_top.spines['left'].set_visible(False)
        ax_right.spines['top'].set_visible(False)
        ax_right.spines['right'].set_visible(False)
        ax_right.spines['bottom'].set_visible(False)
        
        # Remove all ticks for marginal plots
        ax_top.set_xticks([])
        ax_top.set_yticks([])
        ax_right.set_xticks([])
        ax_right.set_yticks([])
        
        # Set background color to match main plot
        ax_top.set_facecolor('white')
        ax_right.set_facecolor('white')
    
    # Add error annotations if requested
    if Err_comment:
        error_file = f'{plot_type}_error.txt' if plot_type in ['energy', 'force'] else 'error.txt'
        err_data = None
        if os.path.exists(error_file):
            try:
                err_data = np.loadtxt(error_file)
                logging.info(f'Loading {plot_type} error annotation data')
            except Exception as e:
                logging.error(f'Failed to load {error_file}: {e}')
        else:
            logging.warning(f'{error_file} file does not exist, error annotations will not be added')
            
        if err_data is not None:
            if err_data.ndim == 1:
                err_data = err_data.reshape(1, -1)
                
            if err_data.shape[1] == 3 and plot_type == 'energy':
                logging.debug('Adding energy error annotations')
                err_img = err_data[:, 0]
                err_point = err_data[:, 1:3]
                # Reduce annotation density to avoid overlap
                max_annotations = 20  # Limit number of annotations
                if len(err_data) > max_annotations:
                    # Sample annotations to avoid overcrowding
                    indices = np.linspace(0, len(err_data)-1, max_annotations, dtype=int)
                    err_img = err_img[indices]
                    err_point = err_point[indices]
                
                for i in range(len(err_img)):
                    ax_main.annotate(f'{int(err_img[i])}', err_point[i], 
                                   textcoords="offset points", xytext=(3,3), 
                                   ha='left', va='bottom', color='red', 
                                   fontsize=max(5, fontsize-8), fontweight='bold',
                                   bbox=dict(boxstyle='round,pad=0.2', fc='white', alpha=0.7))
                                   
            elif err_data.shape[1] == 4 and plot_type == 'force':
                logging.debug('Adding force error annotations')
                err_img = err_data[:, 0]
                err_img_atom = err_data[:, 1]
                err_point = err_data[:, 2:4]
                # Reduce annotation density for force plots
                max_annotations = 15  # Even fewer for force plots due to more data
                if len(err_data) > max_annotations:
                    indices = np.linspace(0, len(err_data)-1, max_annotations, dtype=int)
                    err_img = err_img[indices]
                    err_img_atom = err_img_atom[indices]
                    err_point = err_point[indices]
                
                for i in range(len(err_img)):
                    ax_main.annotate(f'{int(err_img[i])}-{int(err_img_atom[i])}', err_point[i], 
                                   textcoords="offset points", xytext=(2,2), 
                                   ha='left', va='bottom', color='red', 
                                   fontsize=max(4, fontsize-9), fontweight='bold',
                                   bbox=dict(boxstyle='round,pad=0.1', fc='white', alpha=0.6))
    
    try:
        # Keep simple file names as requested
        output_file = f'{plot_type}_error_analysis.png'
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        logging.info(f'Plot saved as {output_file}')
    except Exception as e:
        logging.error(f'Error saving plot: {e}')
    finally:
        plt.close()


def plot_scatter(x, 
                 y, 
                 fontsize, 
                 Err_commennt=True, 
                 colors=['blue','green','orange'], 
):
    """Legacy plot function - kept for compatibility"""
    logging.warning('plot_scatter is deprecated, use plot_scatter_with_marginals instead')
    
    # Determine plot type based on data shape
    plot_type = 'energy' if x.shape[1] == 1 else 'force'
    plot_scatter_with_marginals(x, y, fontsize, plot_type, Err_commennt, colors)

def main(trajname, 
         apps, 
         resource,
         fontsize, 
         data=None, 
         pot=None, 
         er_bar=1.5, 
         ra='Au', 
         cut_img=True, 
         comment=False,
         show_marginals=True
):
    """
    Main function for ML potential error analysis
    
    Parameters:
    -----------
    trajname : str
        Trajectory file path
    apps : str  
        ML software name ('nep' or 'n2p2')
    resource : str
        Data source type ('software' or 'images')
    fontsize : int
        Plot font size
    data : str, dict, or None, optional
        Training results file path(s). Can be:
        - str: single file path (legacy format, auto-detect energy/force)
        - dict: {'energy': 'energy_file.txt', 'force': 'force_file.txt'}
        - None: required when resource='images'
    pot : str, optional
        ML potential file path (required when resource='images')
    er_bar : float, optional
        Error threshold multiplier (default: 1.5)
    ra : str or dict, optional
        Replacement atom type (default: 'Au')
    cut_img : bool, optional
        Whether to save acceptable structures (default: True)
    comment : bool, optional
        Whether to add error annotations (default: False)
    show_marginals : bool, optional
        Whether to show marginal density distributions (default: True)
    """
    try:
        traj = read(trajname, ':')
        logging.info(f'Successfully loaded {len(traj)} structures from {trajname}')
    except Exception as e:
        logging.error(f'Failed to load trajectory file: {e}')
        return
        
    logging.info(f'Starting ML error analysis:')
    logging.info(f'  - Software: {apps}')
    logging.info(f'  - Data source: {resource}')
    logging.info(f'  - Error threshold: {er_bar} × RMSE')
    logging.info(f'  - Analyzing both energy and force data')
    
    try:
        data_normalization(apps, traj, resource, data, pot)
        err_structure_finding(error_bar=er_bar, images=traj, fontsize=fontsize, replace_atom=ra, Cutimg=cut_img, comment=comment, show_marginals=show_marginals)
        logging.info('ML error analysis completed')
    except Exception as e:
        logging.error(f'Error during analysis: {e}')
        raise