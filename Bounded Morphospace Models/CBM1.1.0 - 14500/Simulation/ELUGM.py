from cc3d import CompuCellSetup

from ELUGMSteppables import ELUGMSteppable
CompuCellSetup.register_steppable(steppable = ELUGMSteppable(_frequency = 1))

from ELUGMSteppables import MitosisSteppable
CompuCellSetup.register_steppable(steppable = MitosisSteppable(_frequency = 1))

CompuCellSetup.run()
        
        