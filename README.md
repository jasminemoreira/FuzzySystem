# FuzzySystem

This set of scripts makes up the implementation of a fuzzy system. They can be used to create experiments or production systems. Here is a basic example of their use:

```
# Load libraries
from FuzzyPlots import plot, plotMesh
from FuzzyBasicSets import TriangularSet, SingletonSet, Universe
from FuzzyRelations import FuzzyRelation
from FuzzySystem import FuzzyPartition, FuzzySystem

# Define universe ranges
UTxProdPlan = Universe(-3, 3, term='UTxProdPlan')
UMSE = Universe(0, 0.8, term='MSE')
UACT = Universe(0, 10, term='ACT')

# Define antecedents
TxProdPlan = FuzzyPartition([TriangularSet(-1E10, -2, -1, 1,'GD'),
                    TriangularSet(-2, -1, 0, 1,'MD'), 
                    TriangularSet(-1, 0, 1, 1,'OK'),
                    TriangularSet(0, 1, 2, 1,'ME'),
                    TriangularSet(1, 2, 1E10, 1,'GE')],
                    'TxProdPlan')
plot(UTxProdPlan, TxProdPlan)

MSE = FuzzyPartition([TriangularSet(0, 0, .2, 1,'PE'),
                    TriangularSet(0, .2, .4, 1,'ME'), 
                    TriangularSet(.2, .4, 1E10, 1,'GE')],
                    'MSE')
plot(UMSE, MSE)

# Define consequent
ACT = FuzzyPartition([TriangularSet(0, 0, 5, 1,'BA'),
                    TriangularSet(0, 5, 10, 1,'MA'), 
                    TriangularSet(5, 10, 1E10, 1,'GA')],
                    'Attention')
plot(UACT, ACT)

# Define the base of rules
BR  = [[TxProdPlan.GD, MSE.PE, ACT.GA],
       [TxProdPlan.MD, MSE.PE, ACT.MA],
       [TxProdPlan.OK, MSE.PE, ACT.BA],
       [TxProdPlan.ME, MSE.PE, ACT.MA],
       [TxProdPlan.GE, MSE.PE, ACT.GA],
       [TxProdPlan.GD, MSE.ME, ACT.GA],
       [TxProdPlan.MD, MSE.ME, ACT.GA],
       [TxProdPlan.OK, MSE.ME, ACT.BA],
       [TxProdPlan.ME, MSE.ME, ACT.GA],
       [TxProdPlan.GE, MSE.ME, ACT.GA],
       [TxProdPlan.GD, MSE.GE, ACT.GA],
       [TxProdPlan.MD, MSE.GE, ACT.GA],
       [TxProdPlan.OK, MSE.GE, ACT.GA],
       [TxProdPlan.ME, MSE.GE, ACT.GA],
       [TxProdPlan.GE, MSE.GE, ACT.GA]]

# Create the Fuzzy System
FS = FuzzySystem([UTxProdPlan, UMSE, UACT], BR)

# Define inputs for inference
TxProdPlanInp = SingletonSet(-0.9)
MSEInp = SingletonSet(0.19)
INP = [TxProdPlanInp, MSEInp]

# Map the inputs in their respective universes
UTxProdPlan.mapSet(TxProdPlanInp)
UMSE.mapSet(MSEInp)

# Do the inference
FS.inference(INP, aggregator=9, semantics="Mamdani", defuzzyfication="MoM", precision_mode=False)
print("Conclusão defuzzificada (MoM): ", FS.MoM)
print("Conclusão defuzzificada (Centroid): ", FS.calcCentroid())

# Plot results
plot(UACT, FS.ICONC, show_overlaid=True) 
plot(UACT, FS.getConclusionAsPartition())
```
