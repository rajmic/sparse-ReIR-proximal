# sparse-ReIR-proximal
This repository provides Matlab codes related to the paper "Fast reconstruction of sparse relative impulse responses via second-order cone programming" (authors Rajmic, Koldovský, Daňková, published at WASPAA workshop 2017)

In folder SYNTHETIC, you find Matlab codes showing the algorithms on simple randomly generated data. Just run the file demo_sparse_impulse_response_recovery.m.
In folder SPEECH, the algorithms are applied to the signal described in the paper. Run the file RunME.m or call directly the function experimentSOCP. Please follow the instructions around line 30 to load the sounds. Note that the experiment is performed in 100 trial by default, which may take hours to finish.

In case of trouble please contact us at rajmic@feec.vutbr.cz or zbynek.koldovsky@tul.cz.
