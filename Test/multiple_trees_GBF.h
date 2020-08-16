#pragma once
#ifndef MULTIPLE_TREES_GBF_H
#define MULTIPLE_TREES_GBF_H

#include "kd_tree.h"
#include "ball_tree.h"
#include "GBF.h"

void update_models(vector<model>& our_model_vec, model& our_model);
void build_multiple_trees(vector<Tree>& multi_trees, model& our_model);
void multiple_GBF_iter(int q_id, vector<Tree>& multi_trees, model& our_model);

#endif