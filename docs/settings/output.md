Output format
=====

The base output file is TSV (Tab-Separated Values) format, including the following contents:

| column id              | meaning                                                      |
| ---------------------- | ------------------------------------------------------------ |
| Bac1                   | ID of the 1st bacterium                                      |
| Bac2                   | ID of the 2nd bacterium                                      |
| Growth1                | Growth rate of the 1st bacterium                             |
| Growth2                | Growth rate of the 2nd bacterium                             |
| Intermediate           | The cross-feeding metabolite we used to build the community  |
| EX_{substrate}_1       | Intake rate of given substrate by the 1st bacterium          |
| EX_{substrate}_2       | Intake rate of given substrate by the 2nd bacterium          |
| Glucose_absorption_1   | Intake rate of glucose by the 1st bacterium                  |
| Glucose_absorption_2   | Intake rate of glucose by the 2nd bacterium                  |
| EX_{product}_1         | Production rate of given product by the 1st bacterium        |
| EX_{product}_2         | Production rate of given product by the 2nd bacterium        |
| Cross_feeding_forward  | Cross feeding metabolites from the 1st bacterium to the 2nd bacterium, including their exchange rate |
| Cross_feeding_reverse  | Cross feeding metabolites from the 2nd bacterium to the 1st bacterium, including their exchange rate |
| Bac1_mono_growth       | Growth rate of the 1st bacterium in mono-culture             |
| Bac2_mono_growth       | Growth rate of the 2nd bacterium in mono-culture             |
| DeepCooc_Co_occurrence | Prediction of co-occurrence by DeepCooc                      |
| Interaction_type       | Interaction type of the microbial community                  |
| Total_production       | Cumulative production rate of the microbial community        |
