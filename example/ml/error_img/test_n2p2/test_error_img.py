from qctools.ml.error_img import main

# Updated example without option parameter - now analyzes both energy and force
# Test with marginals and comments disabled
main(trajname='input_data.traj',
     apps='n2p2',
     resource='software',
     fontsize=12,
     data={'energy':'trainpoints.000040.out', 'force':'trainforces.000040.out'},
     pot=None,
     er_bar=2.5,
     ra={'Pt':'Pd', 'H':'He', 'O':'F'},
     cut_img=True,
     comment=False,      # Disable error annotations to reduce clutter
     show_marginals=False  # Enable beautiful KDE marginal distributions
)