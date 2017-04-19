#!/bin/bash

temp=$(grep 'Temp' input_parameters.txt | cut -d ' ' -f 2)
msc=$(grep 'MSC' input_parameters.txt | cut -d ' ' -f 2)
max_bb=$(grep 'max_bb' input_parameters.txt | cut -d ' ' -f 2)
max_cm=$(grep 'max_cm' input_parameters.txt | cut -d ' ' -f 2)
bbrat=$(grep 'bbrat' input_parameters.txt | cut -d ' ' -f 2)
n_meas=$(grep 'N_meas' input_parameters.txt | cut -d ' ' -f 2)
n_a=$(grep 'N_A' input_parameters.txt | cut -d ' ' -f 2)
n_b=$(grep 'N_B' input_parameters.txt | cut -d ' ' -f 2)
l=$(grep 'L' input_parameters.txt | cut -d ' ' -f 2)
n_beads=$(grep 'N_beads' input_parameters.txt | cut -d ' ' -f 2)
g_a=$(grep 'g_A=' input_parameters.txt | cut -d ' ' -f 2)
g_b=$(grep 'g_B' input_parameters.txt | cut -d ' ' -f 2)
g_ab=$(grep 'g_AB' input_parameters.txt | cut -d ' ' -f 2)
m_a=$(grep 'M_A' input_parameters.txt | cut -d ' ' -f 2)
m_b=$(grep 'M_B' input_parameters.txt | cut -d ' ' -f 2)

echo $temp $msc $max_bb $max_cm $bbrat $n_meas > pathintegral_input.txt
echo $n_a $n_b $m_a $m_b $g_a $g_ab $g_ab $l $n_beads > initial_conf_interaction.txt
