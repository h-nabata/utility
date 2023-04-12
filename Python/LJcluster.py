import sys
import os
import numpy as np

epsilon = 1.0;  sigma = 1.7                                          # the parameter of the Lennard-Jones potential
atomnum = 21                                                         # the number of atoms in the system
maxitration = 30;  meta_maxitration = 200                            # the number of iteration in geometry optimization
stepsize = 0.1;  criteria = 1e-4;  stepsize_ene_criteria = 1e-4      # parameters of optimization (@ steepest descent method)
attenuation_mode = 3;  attenuation_scale = 0.95                      # control if the stepsize attenuate in the optimization by the steepest descent method (1 ~ 3)
debug_mode = 0                                                       # control if the value of energy, gradient vectors and stepsize should be written in each iteration
fileoutput_mode = 0                                                  # control if the value of energy, gradient vectors and stepsize should be written to another output file in each run (1 or 2)
dummy_atom_name = "Ar"                                               # specify an element name (dummy) for xyz output
EQlist_mode = 1;  max_explore_runnum = 10                            # control if output of the obtained EQs, and the number of times of exploration trials
outputfilepath = 'C:/Users/'+os.getlogin()+'/Downloads/'             # specify the PATH for output files (a "slash" should be included at the end of the string)


def LJpot_ene(r, atomnum):
    # return L-J potential energy matrix between all atom pairs (atomnum x atomnum)
    with np.errstate(invalid='ignore'):  # ignore the error message "RuntimeWarning: invalid value encountered in subtract"
        a = 4.0 * epsilon * ((sigma/r)**12-(sigma/r)**6)
    return np.nan_to_num(a, nan=0.0).reshape(atomnum, atomnum)  # substitute "nan" to zero value

def LJpot_grad(r, atomnum):
    # return L-J potential gradient matrix between all atom pairs (atomnum x atomnum)
    with np.errstate(invalid='ignore'):  # ignore the error message "RuntimeWarning: invalid value encountered in subtract"
        a = -4.0 * epsilon * (12 * (sigma**12/r**13) - 6 * (sigma**6/r**7))
    return np.nan_to_num(a, nan=0.0).reshape(atomnum, atomnum)  # substitute "nan" to zero value

def norm_check(matrix):
    # return norm of the matrix or vectors  ## https://numpy.org/doc/stable/reference/generated/numpy.linalg.norm.html
    return np.linalg.norm(matrix, ord=2)

# ignore zero division errors
np.seterr(divide='ignore')

# generate initial coordinates randomly (4 x 3 list)
random_range = atomnum ** (1/2)
init_coord_list = random_range * np.random.rand(atomnum, 3) - random_range / 2

# reserve the value of the original stepsize
original_stepsize = stepsize
foundnum = 0
if EQlist_mode >= 1:
    # generate an EQ_list file
    with open(outputfilepath+'LJopt_EQ_list.xyz', mode='w', encoding='utf-8') as fout:
        pass
else:
    max_explore_runnum = 1

# begin optimization
for explore_runnum in range(max_explore_runnum):

    if fileoutput_mode >= 1:
        # generate an output file to record the energy and geometry at the end of each run
        with open(outputfilepath+'LJopt_RUN.xyz', mode='w', encoding='utf-8') as fout:
            pass

    opt_done = 0
    itr_runnum = 0
    metaitr_runnum = 0
    if EQlist_mode >= 1:
        random_range = atomnum ** (1/2)
        init_coord_list = random_range * np.random.rand(atomnum, 3) - random_range * 0.5
    atomnum = len(init_coord_list)
    coord_list = np.array(init_coord_list).reshape(atomnum, 3)

    stepsize = original_stepsize
    tmpene1 = 1; tmpene2 = 1; tmpene3 = 1
    tmpnormval1 = 10;   tmpnormval2 = 10;   tmpnormval3 = 10

    for metaitrnum in range(meta_maxitration):

        if fileoutput_mode >= 2:
            # generate an output file to record the energy and geometry at the end of each meta-iteration
            with open(outputfilepath+'LJopt'+str(metaitrnum)+'.xyz', mode='w', encoding='utf-8') as fout:
                pass

        for itrnum in range(maxitration):
            if debug_mode >= 3:
                print("\n# ITR.", str(itrnum), "coord_list\n", coord_list)

            # calculate distances between all atom pairs (atomnum x atomnum)
            distance_list = np.array([])
            for i in range(atomnum):
                tmp_distance_list = np.array([])
                for j in range(atomnum):
                    tmp_distance = 0.0
                    for k in range(3):
                        tmp_distance = tmp_distance + (coord_list[i][k] - coord_list[j][k]) ** 2
                    tmp_distance_list = np.append(tmp_distance_list, np.sqrt(tmp_distance))
                distance_list = np.append(distance_list, tmp_distance_list, axis=0)
            distance_list = distance_list.reshape(atomnum, atomnum)
            ene_matrix = LJpot_ene(distance_list, atomnum)
            grad_matrix = LJpot_grad(distance_list, atomnum)
            if debug_mode >= 3:
                print("distance_list\n", distance_list.reshape(atomnum, atomnum))
                print("ene_matrix\n", ene_matrix)
                print("grad_matrix\n", grad_matrix)

            # difine xyz vectors between all atom pairs (atomnum x atomnum)
            tmp_vector_list = np.array([])
            for i in range(atomnum):
                for j in range(atomnum):
                    tmp_vector = np.array([])
                    for k in range(3):
                        tmp_vector = np.append(tmp_vector, (coord_list[i][k] - coord_list[j][k]))
                    tmp_vector_list = np.append(tmp_vector_list, tmp_vector, axis=0)
            vector_list = tmp_vector_list.reshape(atomnum, atomnum, 3)
            if debug_mode >= 3:
                print("vector_list\n", vector_list)

            # calculate force vectors (atomnum x 3)
            force_list = np.array([])
            for i in range(atomnum):
                tmp_force_list = np.array([0.0, 0.0, 0.0])
                for j in range(atomnum):
                    if distance_list[i][j] < 1e-1:
                        pass
                    else:
                        tmp_force_list = tmp_force_list + vector_list[i][j] * grad_matrix[i][j] / distance_list[i][j]
                force_list = np.append(force_list, tmp_force_list)
            force_list = force_list.reshape(atomnum, 3)
            normval = norm_check(force_list)
            if debug_mode >= 3:
                print("force_list\n", force_list)
                print("norm :", normval)

            # calculate the energy
            optene = 0.0
            for i in range(atomnum):
                for j in range(atomnum):
                    optene = optene + ene_matrix[i][j]
            tmpene1 = tmpene2;   tmpene2 = tmpene3;   tmpene3 = optene

            # output section (option)
            if fileoutput_mode >= 1:
                if fileoutput_mode >= 2:
                    with open(outputfilepath+'LJopt'+str(metaitrnum)+'.xyz', mode='a', encoding='utf-8') as fout:
                        fout.write(str(atomnum)+'\n')
                        fout.write('# ITR. '+str(itrnum)+', E = '+str(optene/2)+'\n')
                        for i in range(len(coord_list)):
                            line = dummy_atom_name + " " + str(coord_list[i][0]) + " " + str(coord_list[i][1]) + " " + str(coord_list[i][2])
                            fout.write(line+'\n')
                        fout.write('\n')

                if itrnum == maxitration - 1:
                    with open(outputfilepath+'LJopt_RUN.xyz', mode='a', encoding='utf-8') as fout2:
                        fout2.write('# RUN. '+str(metaitrnum)+', E = '+str(optene/2)+'\n')
                        for i in range(len(coord_list)):
                            line = dummy_atom_name + " " + str(coord_list[i][0]) + " " + str(coord_list[i][1]) + " " + str(coord_list[i][2])
                            fout2.write(line+'\n')
                        fout2.write('\n')

            # cahnge the stepsize (option)
            if attenuation_mode == 1:
                # cause a attenuation in stepsize (type 1)
                stepsize = stepsize * (maxitration - itrnum) / maxitration
            elif attenuation_mode >= 2:
                # cause a attenuation in stepsize (type 2)
                if tmpene1 + tmpene2 + tmpene3 < stepsize_ene_criteria:
                    if debug_mode >= 3:
                        print("The criteria met.\nstepsize :", stepsize, "->", stepsize * attenuation_scale)
                        print("stepsize_ene_criteria :", stepsize_ene_criteria, "->", stepsize_ene_criteria * attenuation_scale)
                    stepsize = stepsize * attenuation_scale
                    stepsize_ene_criteria = stepsize_ene_criteria * attenuation_scale

            if normval < criteria:
                opt_done = 1
                itr_runnum = itrnum
                break

            # update the coordinates
            coord_list = coord_list - (stepsize / normval) * force_list

        if attenuation_mode >= 3:
            # cause a jump (recovery) in stepsize
            pre_normval = tmpnormval1 + tmpnormval2 + tmpnormval3
            tmpnormval1 = tmpnormval2;  tmpnormval2 = tmpnormval3;  tmpnormval3 = normval
            sub_normval = tmpnormval1 + tmpnormval2 + tmpnormval3
            if (pre_normval - sub_normval) / pre_normval < 1e-4 :
                stepsize = original_stepsize
        if debug_mode >= 1:
            print("norm, stepsize :", normval, stepsize)
            if debug_mode >= 2:
                print("force_list\n", force_list)

        if opt_done == 1:
            metaitr_runnum = metaitrnum
            break

    # if an optimized structure is found
    if opt_done == 1:
        print("\n--------------------------\ntrial", str(explore_runnum), "; Optimized geometry ( RUN :", str(metaitr_runnum), "/ ITR :", str(itr_runnum),")\n", coord_list)
        # print("ene_matrix\n", ene_matrix)
        optene = 0.0
        for i in range(atomnum):
            for j in range(atomnum):
                optene = optene + ene_matrix[i][j]
        print("E = ", optene/2)
        print("Force vectors\n", force_list)
        if fileoutput_mode >= 1:
            with open(outputfilepath+'LJopt_RUN.xyz', mode='a', encoding='utf-8') as fout2:
                fout2.write('\n--------------------------\nOptimized geometry ( RUN :'+str(metaitrnum)+" / ITR :"+str(itr_runnum)+'), E = '+str(optene/2)+'\n')
                for i in range(len(coord_list)):
                    line = dummy_atom_name + " " + str(coord_list[i][0]) + " " + str(coord_list[i][1]) + " " + str(coord_list[i][2])
                    fout2.write(line+'\n')
                fout2.write('\n')
        if EQlist_mode >= 1:
            with open(outputfilepath+'LJopt_EQ_list.xyz', mode='a', encoding='utf-8') as fout3:
                fout3.write(str(atomnum)+'\nEQ '+str(foundnum)+' ; E = '+str(optene/2)+'\n')
                for i in range(len(coord_list)):
                    line = dummy_atom_name + " " + str(coord_list[i][0]) + " " + str(coord_list[i][1]) + " " + str(coord_list[i][2])
                    fout3.write(line+'\n')
                fout3.write('\n')
        foundnum += 1
    else:
        print("\n--------------------------\nOptimization did not completed\n")
        print("Force vectors\n", force_list)
