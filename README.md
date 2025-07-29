# Post-selection Inference for Conformal Prediction  
**[arXiv:2304.06158](https://arxiv.org/abs/2304.06158)**  

## Overview

Conformal prediction is a powerful tool for providing **distribution-free, finite-sample valid prediction sets** for black-box machine learning models. However, traditional conformal inference requires a **pre-specified, data-independent miscoverage level**.

In many practical situations, analysts may want to **adapt the miscoverage level post hoc**—after observing the prediction sets—to trade off **coverage for precision**.

## Motivation

For example, in binary classification:
- A conformal prediction set with 95% coverage might frequently include **both classes**, making it uninformative.
- The analyst may wish to examine **80% coverage** instead, to obtain more precise (i.e., smaller) sets that still retain meaningful uncertainty quantification.

Such **data-driven selection of coverage levels** introduces a **post-selection inference problem**, where standard conformal guarantees no longer apply.

## Contribution

We develop a framework for **simultaneous conformal inference**, enabling:
- **Valid post-selection inference**: Simultaneous coverage guarantees hold for **all** miscoverage levels in a pre-specified range.
- **Finite-sample guarantees** under the i.i.d. assumption.
- **User flexibility**: Practitioners can adaptively choose the desired coverage level *after* seeing the prediction sets—by any criterion (e.g., size or informativeness)—without compromising statistical validity.

## Key Features

- 📏 **Simultaneous coverage guarantee**:  
  For a range of miscoverage levels \(\alpha \in I\), we ensure  
  \[
  \Pr\left(\text{prediction set covers } Y_i, \ \forall \alpha \in I \right) \geq 1 - \delta.
  \]

- 🔄 **Post-hoc flexibility**:  
  Analysts can trade off **prediction set size** and **coverage** without retraining or violating inference guarantees.

- 📊 **Practical utility**:  
  Useful in classification, regression, and general black-box ML tasks where **prediction set quality matters**.

## Paper

See the full paper here:  
📄 [arXiv:2304.06158](https://arxiv.org/abs/2304.06158)

## Citation

```bibtex
@article{sarkar2023postselection,
  title={Post-selection Inference for Conformal Prediction: Trading off Coverage for Precision},
  author={Siddhaarth Sarkar and Arun Kumar Kuchibhotla},
  journal={arXiv preprint arXiv:2304.06158},
  year={2023}
}
