import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import re

plt.rcParams.update({
    "font.family": "STIXGeneral",
    "mathtext.fontset": "stix",
    "font.size": 16
})

def parse_filband(feig, npl=10):
    with open(feig, 'r') as f:
        lines = f.readlines()

    header = lines[0].strip()
    shape = re.split('[,=/]', header)
    nbnd = int(shape[1])
    nks = int(shape[3])
    eig = np.zeros((nks, nbnd), dtype=np.float32)

    div = nbnd // npl + (0 if nbnd % npl == 0 else 1)
    kinfo = []
    for index, value in enumerate(lines[1:]):
        value = value.strip()
        quotient = index // (div + 1)
        remainder = index % (div + 1)
        if remainder == 0:
            kinfo.append(value)
        else:
            vals = re.split(r'\s+', value)
            a = (remainder - 1) * npl
            b = a + len(vals)
            eig[quotient][a:b] = vals

    return eig, nbnd, nks, kinfo

def read_pdos(file_name):
    data = np.loadtxt(file_name, comments='#')
    energies = data[:, 0]
    pdos = data[:, 1]
    return energies, pdos

def read_pdos_tot(file_name):
    data = np.loadtxt(file_name, comments='#')
    energies_tot = data[:, 0]
    pdos_tot = data[:, 2]
    return energies_tot, pdos_tot

def draw_proj_band(proj_file, bd_file, fig_file):
    eig, nbnd, nks, kinfo = parse_filband(bd_file)

    E_Fermi = -0.2217  # Fermi level
    ymin, ymax = -5, 5
    lw = 0.5

    # Orbital groupings (update as needed)
    oo = [
    [0,1,2,3,4,5,6,7,8,9],  # Nb 
    [10,11,12,13,14,15,16,17,18,19],  # Mo 
    [20,21,22,23],  # C 
    [24,25,26,27,29,30,31,28],  # O 
]

    colors = ['#002850', '#e8b950', '#bbbbbb', '#d41d4e']
    labels = ['Nb', 'Mo', 'C', 'O']
    scale = 200

    x_ticks = [0,  100, 200, 300, 400]
    x_labels = [ 'K', r'$\Gamma$', 'M', 'K', r'$\Gamma$' ]

    # Prepare PDOS
    energies_Nb, pdos_Nb = read_pdos('atom_Nb.dat')
    energies_Mo, pdos_Mo = read_pdos('atom_Mo.dat')
    energies_C, pdos_C = read_pdos('atom_C.dat')
    energies_O, pdos_O = read_pdos('atom_O.dat')
    energies_tot, pdos_tot = read_pdos('NbMoCO2.k.pdos_tot')

    # Align energies to Fermi level
    energies_Nb -= E_Fermi
    energies_Mo -= E_Fermi
    energies_C -= E_Fermi
    energies_O -= E_Fermi
    energies_tot -= E_Fermi
    eig -= E_Fermi

    # Read projected bands
    with open(proj_file, 'r') as filproj:
        for _ in range(16):
            filproj.readline()

        nlorb, nks_read, nbnd_read = map(int, filproj.readline().split()[:3])
        pjmat = np.zeros((nlorb, nks, nbnd), dtype=np.float32)

        filproj.readline()
        for i in range(nlorb):
            filproj.readline()
            for j in range(nks):
                for k in range(nbnd):
                    pjmat[i, j, k] = float(filproj.readline().split()[2])

    # Plotting: BANDS + PDOS
    fig = plt.figure(figsize=(6, 5))
    gs = gridspec.GridSpec(1, 2, width_ratios=[3, 1], wspace=0.05)

    ax = plt.subplot(gs[0])
    for i in range(nbnd):
        ax.plot(np.arange(nks), eig[:, i], color='gray', linewidth=lw)
    ax.axhline(0, color='gray', linestyle='--', linewidth=1)
    for vline in x_ticks:
        ax.axvline(x=vline, ymin=0, ymax=1, linewidth=lw, color='black')

    # Projections
    for i, orb in enumerate(oo):
        ax.scatter(-1, ymin-1, 20, c=colors[i], alpha=0.8, label=labels[i], marker='.', edgecolors='none')
        for k in range(nbnd):
            s_of_o = np.sum(pjmat[orb, :, k], axis=0)
            ax.scatter(np.arange(nks), eig[:, k], s=scale * s_of_o**2, c=colors[i], alpha=0.5, marker='.', edgecolors='none')

    ax.set_xlim([0, nks - 1])
    ax.set_ylim([ymin, ymax])
    ax.set_ylabel('E - E$_F$ (eV)', fontsize=16, labelpad=2)
    ax.set_xticks(x_ticks)
    ax.set_xticklabels(x_labels)


    # Plot PDOS
    ax_pdos = plt.subplot(gs[1], sharey=ax)
    ax_pdos.plot(pdos_Nb, energies_Nb, label='Nb', color=colors[0])
    ax_pdos.plot(pdos_Mo, energies_Mo, label='Mo', color=colors[1])
    ax_pdos.plot(pdos_C, energies_C, label='C', color=colors[2])
    ax_pdos.plot(pdos_O, energies_O, label='O', color=colors[3])
    ax_pdos.plot(pdos_tot, energies_tot, label='Tot', color='black', linestyle='--')

    ax_pdos.axhline(0, color='gray', linestyle='--', linewidth=1)
    ax_pdos.set_xlim([0, 7])
    ax_pdos.set_xticks([])
    ax_pdos.set_xlabel('PDOS', fontsize=14)
    ax_pdos.tick_params(axis='y', labelleft=False)
    
    ax.legend(
        scatterpoints=1, markerscale=2.0, loc='upper right',
        handlelength=1.2, handletextpad=0.4, borderpad=0.3,
        labelspacing=0.3, fontsize=10, frameon=True, fancybox=True
    )
    ax_pdos.legend(
    loc='upper right',
    fontsize=10,
    frameon=True,
    fancybox=True,
    borderpad=0.3,
    labelspacing=0.3,
    handlelength=1.2,
    handletextpad=0.4
)
    #fig.text(0.05, 0.92, '(a)', fontsize=16, va='top')
    #fig.text(0.72, 0.92, '(b)', fontsize=16, va='top')
    plt.savefig(fig_file, dpi=300, bbox_inches='tight')
    plt.show()

if __name__ == '__main__':
    draw_proj_band("NbMoCO2proj.out.projwfc_up", "NbMoCO2.band", "Proj_Band_PDOS.pdf")
