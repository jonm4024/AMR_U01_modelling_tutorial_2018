# test = "5000_vs_50000_samples"  # current test name
# mcmc_reduced_dict_p1[test] contains the data frame that includes the data from sampling.
# For instance, this is a comparison between 5000 and 50000, therefore
# mcmc_reduced_dict_p1[test] contains 5000 for the varying parameter.

#mcmc_reduced_dict_p2[test]
# mcmc_reduced_dict_p2[test] always contains the data frame that includes the
# data from sampling at 50000 for the varying parameter.


# Code necessary to enable calculate_pvalue_permutation to execute without crashing.
# Takes the intersection of reaction df's to calculate stats on. Otherwise
# calculate_pvalue_permutation crashes, complaining that for "z = cond1 - cond2",
# "shapes" are not the same size.
# set1=set(mcmc_reduced_dict_p1[test].columns.tolist())
# set2=set(mcmc_reduced_dict_p2[test].columns.tolist())
# intersect=list(set1&set2)
# intersect_dict_p1 = mcmc_reduced_dict_p1[test][intersect]
# intersect_dict_p2 = mcmc_reduced_dict_p2[test][intersect]
# flux_stats_dict_p1_p2 = {}

# for test in mcmc_reduced_dict_p1.keys():
        # print "test:", test

        # flux_stats_dict_p1_p2[test],
        # delta_mean_df, delta_median_df,
        # delta_std_df = get_mean_std_pvalue(intersect_dict_p1,
                                            # intersect_dict_p2, m, 1, 2)
                                            
import numpy as np
import pandas as pd

def calculate_pvalue_permutation(data_1_I, data_2_I, n_permutations_I = 10, n_resamples_I = 10):
    '''calculate the pvalue of two data by determining
    the lack of overlap between sample points using a permutation test.
    If the sample points of the data sets is not equal,
    a subset of samples of matching length is resampled from the larger data set'''

    data_1 = None; #sample set with fewer points
    data_2 = None; #sample set with more points
    n_resamples = 0;



    # check the length of data_1 and data_2
    if len(data_1_I)>len(data_2_I):
        data_1=np.array(data_2_I);
        data_2=np.array(data_1_I);
        n_resamples=n_resamples_I;
    elif len(data_1_I)<len(data_2_I):
        data_1=np.array(data_1_I);
        data_2=np.array(data_2_I);
        n_resamples=n_resamples_I;
    else:
        data_1=np.array(data_1_I);
        data_2=np.array(data_2_I);

    n_samples_min = len(data_1);

    vals = []
    for i in range(0,n_permutations_I):
        if n_resamples==0:
            cond1 = np.random.permutation(data_1)
            cond2 = np.random.permutation(data_2)
            z = cond1 - cond2
            x = len(z[z>0]) + 1
            y = len(z[z<0]) + 1
            k = min(x,y)
            vals.append(k)
        else:
            cond1 = np.random.permutation(data_1)
            cond2 = np.random.permutation(data_2)
            for resample in range(n_resamples):
                cond2_int = np.random.randint(0,n_samples_min);
                z = cond1 - cond2[cond2_int]
                x = len(z[z>0]) + 1
                y = len(z[z<0]) + 1
                k = min(x,y)
                vals.append(k)
    p = np.mean(vals)/len(data_1)*2
    return p;


def get_mean_std_pvalue(sample1, sample2, m, state1, state2):
    ''' takes in two dataframes that are sampled flux distributions (sample1 and sample2 are two sampled states)
    m is an M model and state1 and state2 are string names of the sampled states. It outputs a new dataframe with
    statistical measures such as means, medians and std's for each state and the pvalue which can be used to determine
    how much overlap there is between two flux distributions of a reaction in each respective state'''

    out = []
    delta_mean_list = []
    delta_median_list = []
    delta_std_list = []

    for r in m.reactions:
        r = r.id
        if r in sample1.columns:
            cond1 = sample1[r]
        else:
            cond1 = pd.Series(np.zeros(10000))
        if r in sample2.columns:
            cond2 = sample2[r]
        else:
            cond2 = pd.Series(np.zeros(10000))


        vals = []
        for i in range(0,10):
            cond1 = cond1.reindex(np.random.permutation(cond1.index))
            cond2 = cond2.reindex(np.random.permutation(cond2.index))
            z = cond1 - cond2
            x = len(z[z>0]) + 1
            y = len(z[z<0]) + 1
            k = min(x,y)
            vals.append(k)
        p = np.mean(vals)/len(cond1)*2  # Previous p value calculation
        # p = calculate_pvalue_permutation(sample1, sample2)  #Douglas' p value calculation.
        c1mean = cond1.mean()
        c2mean = cond2.mean()
        c1median = cond1.median()
        c2median = cond2.median()
        c1std = cond1.std()
        c2std = cond2.std()

        delta_mean = abs(c1mean-c2mean)
        delta_median = abs(c1median-c2median)
        delta_std = abs(c1std-c2std)

        out.append({'reaction':r,'p%dmean'%state1:c1mean, 'p%dmean'%state2:c2mean,'p%dmedian'%state1:c1median, 'p%dmedian'%state2:c2median,'p%dstd'%state1:c1std,'p%dstd'%state2:c2std,'pval':p})
        delta_mean_list.append({'r':r, 'mean_delta_5000':delta_mean})
        delta_median_list.append({'r':r, 'median_delta_5000':delta_median})
        delta_std_list.append({'r':r, 'std_delta_5000':delta_std})

    out = pd.DataFrame(out)
    delta_mean_df = pd.DataFrame(delta_mean_list)
    delta_median_df = pd.DataFrame(delta_median_list)
    delta_std_df = pd.DataFrame(delta_std_list)

    return out, delta_mean_df, delta_median_df, delta_std_df
