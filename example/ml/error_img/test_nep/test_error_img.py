from qctools.ml.error_img import main

# Updated example without option parameter - now analyzes both energy and force
# Test with marginals and comments disabled
main(trajname='train.xyz',
     apps='nep',
     resource='images',
     fontsize=12,
     data=None,
     pot='nep.txt',
     er_bar=2.5,
     ra={'Pt':'Pd', 'H':'He', 'O':'F'},
     cut_img=True,
     comment=False,      # Disable error annotations to reduce clutter
     show_marginals=True  # Enable beautiful KDE marginal distributions
)