# HLS: A High-Performance Local Search Method for Qubit Mapping Problem

*HLS* is a high-performance local search algorithm for solving qubit mapping problem. This repository contains the source code of *HLS*, the coupling graph of NISQ devices and the circuit sets used in the experiment, experiment result files of *HLS* and other competitors and ablation result files.

## Instructions for Building *HLS*

```bash
sh build.sh	
```

By executing the script file, the user can build an executable of *HLS* in the repository root directory. Please note that this script should be run on 64-bit GNU/Linux OS. *HLS* is cross-platform, for other OS, please follow the *Makefile* in the *src* directory to build it.

## Instructions for Running *HLS*

The command to run *HLS* once is as follows:

```bash
./HLS [RANDOM_SEED] [RESULT_FILE_PATH] [DEVICE_NAME] [DELTA] [THETA] [LAMBDA] [QASM_FILE_PATH]
```

DEVICE_NAME is the parameter used to specify the quantum device, and its enumerated values are listed in the following table:
| value      | device in paper |
| ----------- | ----------- |
| `TOKYO`      | IBM Tokyo       |
| `GUADALUPE`      | IBM Guadalupe       |
| `ROCHESTER`   | IBM Rochester        |
| `torino`      | IBM Torino       |
| `SYCAMORE`   | Google Sycamore        |

## Example Command for Running *HLS*

```bash
./HLS 0 ./ TOKYO 5 5 500 ./benchmarks/RevLib/mod5mils_65.qasm
```

Running this command will call the `HLS` program, using `0` as random seed and set the hyperparameters $\beta=5$, $\alpha=5$, $\varphi=500$, map `./benchmarks/RevLib/mod5mils_65.qasm` quantum circuit on the `IBM Tokyo` device, and store the resulting in the `./` directory. For this example, the result file name will be `mod5mils_65.qasm_log.txt`.

## Directories

The directory `device/` includes the coupling graph of NISQ devices.

The directory `src/` includes the implementation of *HLS*. 

The directory `benchmarks/` contains all circuit sets' qasm files. 

The directory `experiment/` contains the result files. 

- compare: provides the original experimental results of *HLS*, *EffectiveQM*, *ILS*, *Qiskit* and *Tket* on all circuit sets and devices.
  - Circuit sets: RevLib, RW, QV and QUEKNO.
  - Devices: TOKYO/tok, GUADALUPE/gua, ROCHESTER/roc, SYCAMORE/syc and torino/tor.
- ablation: provides data for ablation analysis experiments.
  * Ablation versions: Alt1 (bridge-balanced scoring functions ablation), Alt2 (potential-guided tie-braking rule ablation) and Alt3 (iterated initial mapping adjustment ablation).
- analytical: provides analytical experiment results mentioned in Overview section.
  * bridge: the ratio of bridge operation selections to the total number of decisions. EQM*: a version of *EffectiveQM* which only changed the scoring function to $brs$
  * CV: Coefficient of Variation.
  * tie_rate: the proportion of decisions where the highest score is shared by multiple operations. 
- hyperpara: provides hyper-parameter analysis results shown in Appendix.
  * beta: controls the number of iterations in a single search, values are set to 1, 3, 5, 7, and 9.
  * alpha: controls the number of overall search iterations, values are set to 50, 100, 500, and 1000.
  * varphi: controls the depth of scoring, values are set to 1, 3, 5, 7, and 9.

## Main Competitors of *HLS*

As mentioned in the paper, the main competitors of *HLS* are *EffectiveQM*, *ILS*, *Qiskit* and *Tket*. Their codes can be obtained from the following link.

- *EffectiveQM*: https://github.com/chuanluocs/EffectiveQM

- *ILS*: https://github.com/joyofly/ILS-QuantumCircuitMapper

- *Qiskit*: https://github.com/Qiskit/qiskit-terra

- *Tket*: https://github.com/CQCL/tket