# Config for the joint upper-triangle Ising PL pipeline.
# Reuses the same simulated data as the symmetrizing approach (same setting_tag / grid).

source("F:/ising model thesis/-Fourth-Corner-Ising-Model-/Code/ising_config.r")

ising_joint_code_dir = file.path(project_root, "Code", "Ising_joint")
# Estimates / tables live alongside the node-wise Ising outputs, with a _joint tag.
est_tag_joint = paste0(setting_tag, "_joint")

# Adaptive-lasso init for the joint design. "ridge" is safer under separation.
init_method_joint = "unpenalized"
lambda_rule_joint = "lambda.min"
gamma_value_joint = 1
