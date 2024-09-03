Output format
=====

The base output file is TSV (Tab-Separated Values) format, including the following contents:

| column id             | meaning                                                      |
| --------------------- | ------------------------------------------------------------ |
| Bac1                  | id of the 1st bacterium                                      |
| Bac2                  | id of the 2nd bacterium                                      |
| Growth1               | growth rate of the 1st bacterium                             |
| Growth2               | growth rate of the 2nd bacterium                             |
| EX_{substrate}_1      | intake rate of given substrate by the 1st bacterium          |
| EX_{substrate}_2      | intake rate of given substrate by the 2nd bacterium          |
| glucose_absorption_1  | intake rate of glucose by the 1st bacterium                  |
| glucose_absorption_2  | intake rate of glucose by the 2nd bacterium                  |
| EX_{product}_1        | production rate of given product by the 1st bacterium        |
| EX_{product}_2        | production rate of given product by the 2nd bacterium        |
| cross_feeding_forward | cross feeding metabolites from the 1st bacterium to the 2nd bacterium, including their exchange rate |
| cross_feeding_reverse | cross feeding metabolites from the 2nd bacterium to the 1st bacterium, including their exchange rate |
| Bac1_single_growth    | growth rate of the 1st bacterium in mono-culture             |
| Bac2_single_growth    | growth rate of the 2nd bacterium in mono-culture             |
|                       |                                                              |
|                       |                                                              |
