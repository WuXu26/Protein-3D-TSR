import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from collections import Counter
import argparse
import glob
import ntpath
import re

from mpl_toolkits.mplot3d import Axes3D
#import matplotlib as mpl
#mpl.use('Agg')
import matplotlib

__author__ = "Venkata Sarika Kondra"

__version__ = "1.0.1"
__maintainer__ = "Venkata Sarika Kondra"
__email__ = "c00219805@louisiana.edu"


parser = argparse.ArgumentParser(description='Similarities Comparison before and after Size Gap.')
parser.add_argument('--sample_name', '-sample', metavar='sample_name', default='sample_EDD', \
    help='Name of the sample on which this script should be run.')
parser.add_argument('--path', '-path', metavar='path', \
    default='/media/c00219805/Elements1/Research/Hierarchical_Classification/Database/', \
    help='Directory of input sample and other files.')
parser.add_argument('--bin_folder', '-bin_folder', metavar='bin_folder', \
    default='theta29_dist35', \
    help='Directory of input sample and other files.')


def plot_group_analysis(df, ax1, group_name): 
    df = df[['rounded_key_percent']].reset_index()
    df_percent_grouping = df.groupby(['rounded_key_percent']).agg(['count']).reset_index()
    df_percent_grouping.columns = ['key_percent', 'count']
    
    color = 'tab:red'
    ax1.set_xlabel("Percentage of Proteins in only {} group(%)".format(group_name))
    ax1.set_ylabel('Keys Count in {}'.format(group_name), color=color)
    p1,  = ax1.plot(df_percent_grouping['key_percent'].values, 
                    df_percent_grouping['count'].values, color=color, label = group_name)
    ax1.tick_params(axis='y', labelcolor=color)
    #plt.legend(loc = 'best')
    #plt.show()
    return ax1, p1

def plot_unique_keys(df, group_name, t_keys, u_keys, top_n, x_axis_step, y_axis_step):    
    df = df[['rounded_key_percent']].reset_index()
    df_percent_grouping = df.groupby(['rounded_key_percent']).agg(['count']).reset_index()
    df_percent_grouping.columns = ['key_percent', 'count']
    top_n_percents = sorted(list(set(df_percent_grouping['key_percent'].values)))[-top_n:]
    max_percent = np.max(df_percent_grouping['key_percent'].values)
    max_count = np.max(df_percent_grouping['count'].values)
    no_keys = df_percent_grouping.loc[df_percent_grouping['key_percent'] == max_percent, 'count'].values[0]
    plt.bar(df_percent_grouping['key_percent'].values, df_percent_grouping['count'].values)
    plt.xlabel("Percentage of Proteins of only {} group(%)".format(group_name))
    plt.ylabel("Number of Keys")
    plt.title("Distribution of Unique keys in {} group \nw.r.t its Proteins Percentage".\
              format(group_name, group_name))
    plt.ylim(0,max_count + y_axis_step)
    ymin, ymax = plt.gca().get_ylim()
    xmin, xmax = plt.gca().get_xlim()
    x_start = (xmax - xmin)/2
    plt.text(x_start - x_axis_step,  max_count, "#Total Keys: " + r"$\bf{" + str(t_keys) + "}$"  )
    plt.text(x_start- x_axis_step,  (max_count - y_axis_step/2), "#Unique Keys: " + r"$\bf{" + str(u_keys) + "}$" )
    plt.text(x_start- x_axis_step, (max_count - ((3*y_axis_step)/2)), "Max. Percent containing \nunique common keys: " 
             + r"$\bf{" + str(round(np.max(max_percent), 0)) + "}$" + "%")
             
    plt.text(x_start-x_axis_step, (max_count - ((4*y_axis_step)/2)), "#Keys at Max. Percent: "+ r"$\bf{" + str(no_keys) + "}$" )
    for x,y in zip(df_percent_grouping['key_percent'].values, df_percent_grouping['count'].values):
        label = "{}".format(y.round())
        plt.annotate(label, # this is the text
                     (x,y), # this is the point to label
                     textcoords="offset points", # how to position the text
                     xytext=(1,3), # distance from text to points (x,y)
                     ha='center') # horizontal alignment can be left, right or center
    plt.savefig(os.path.join(args.path, args.sample_name, args.bin_folder, 'group_analysis', 
    	"plot_unique_keys_{}.png".format(group_in_test)), bbox_inches='tight')
    return top_n_percents

def generate_two_sided_plot(group_in_test, groups):
    #print(df[df['rounded_key_percent'] == 68.0])
    plots = []
    fig, ax1 = plt.subplots()
    ax1, p1 = plot_group_analysis(df, ax1, groups_display_name[group_in_test])
    plots.append(p1)
    #plt.show()
    ax2 = ax1.twinx()
    ax2.set_ylabel('Average Keys Percentage in Other Groups') 
    percent_values_group = list(set(df['rounded_key_percent'].values))
    all_others_together = {}
    all_files_in_others = 0
    for group in groups:
        if group != group_in_test:
            other_groups_perc = []
            df_g = pd.read_csv(os.path.join(args.path, args.sample_name, args.bin_folder, 
            	'keys_percents_{}_0.csv'.format(group)), header = 0)
            all_files_in_others = all_files_in_others + list(df_g['total_files'].values)[0]
            print(len(df_g))    
            for percent in percent_values_group:
                perc_group_keys = df[df['rounded_key_percent'] == percent]['key'].values
                df_g_perc = df_g[df_g['key'].isin(perc_group_keys)]
                other_groups_perc.append((percent, df_g_perc['key_percent'].mean()))
                if percent in all_others_together.keys():
                    all_others_together[percent] = Counter(all_others_together[percent]) + \
                        Counter(df_g_perc.set_index('key')['#files_present'].to_dict())
                else:
                    all_others_together[percent] = df_g_perc.set_index('key')['#files_present'].to_dict()
                #print(all_others_together, percent)
            df_other_groups = pd.DataFrame(other_groups_perc, columns = ['group_percent', 'other_group_percent_mean'])
            p, = ax2.plot(df_other_groups['group_percent'].values, df_other_groups['other_group_percent_mean'].values, 
                     '--', label = groups_display_name[group])
            plots.append(p)
    ax2.tick_params(axis='y')
    ax1.legend(plots, [l.get_label() for l in plots], loc = 'upper center')
    fig.tight_layout()  # otherwise the right y-label is slightly clipped
    plt.savefig(os.path.join(args.path, args.sample_name, args.bin_folder, 'group_analysis', 
    	"plot_1_{}.png".format(group_in_test)), bbox_inches='tight')
    plt.show()
    plt.close()
    return all_others_together, all_files_in_others

def generate_all_others_together_df(all_others_together, all_files_in_others, group_percent_threshold):
    df_all_others_together = pd.DataFrame(columns = ['key','others_files_count', 'group_percent'])
    for percent, v in all_others_together.items():
        df_f = pd.DataFrame(list(v.items()), columns = ['key','others_files_count'])
        df_f['group_percent'] = percent
        df_all_others_together = df_all_others_together.append(df_f)
    df_all_others_together['all_others_percent'] = 100*df_all_others_together['others_files_count'].div(all_files_in_others)
    #df_all_others_together['diff'] = df_all_others_together['group_percent']-df_all_others_together['all_others_percent']
    df_all_others_together['percent_ratio1'] = df_all_others_together['group_percent']/df_all_others_together['all_others_percent']
    df_all_others_together['percent_ratio1']  = df_all_others_together['percent_ratio1'].apply(lambda x : round(x, 2))
    #df_all_others_together['percent_ratio'] = df_all_others_together['diff'].abs()/df_all_others_together['group_percent']
    #df_all_others_together['percent_ratio']  = df_all_others_together['percent_ratio'].apply(lambda x : round(x, 2))
    df_all_others_together['percent_ratio3']  = df_all_others_together.apply(lambda x : min(x['group_percent'], 
                                                                            x['all_others_percent'])/max(x['group_percent'], 
                                                                                                    x['all_others_percent']), axis =1)
    df_all_others_together['percent_ratio3']  = df_all_others_together['percent_ratio3'].apply(lambda x : round(x, 2))
    df_all_others_together = df_all_others_together.sort_values(['percent_ratio1', 'group_percent'], 
                                                                ascending = False)
    print(df_all_others_together.head(5))
    df_all_others_together = df_all_others_together[\
                                    df_all_others_together['group_percent']> group_percent_threshold]

    #df_all_others_together.to_csv('/media/c00219805/Elements1/Research/Hierarchical_Classification//Database/sample_DD/Class_Analysis/all_others_together_{}_0.csv'.format(group), header = 0)

    plt.hist(df_all_others_together['percent_ratio1'])
    plt.show()
    plt.savefig(os.path.join(args.path, args.sample_name, args.bin_folder, 'group_analysis', 
    	"plot_2_{}.png".format(group_in_test)), bbox_inches='tight')
    plt.close()
    plt.hist(df_all_others_together['percent_ratio3'])
    plt.show()
    plt.savefig(os.path.join(args.path, args.sample_name, args.bin_folder, 'group_analysis', 
    	"plot_3_{}.png".format(group_in_test)), bbox_inches='tight')
    plt.close()
    return df_all_others_together


def generate_allowable_features(df_all_others_together, group_in_test, allowable_feature_ratio):
    # print(sorted(list(set(df_all_others_together['percent_ratio1'].values))))
    print(df_all_others_together[df_all_others_together['percent_ratio1'].isin([allowable_feature_ratio])].head(10))
    print(df_all_others_together[(df_all_others_together['group_percent']> 50 )\
                                 & (df_all_others_together['all_others_percent']< 10 )\
                          ].sort_values(['group_percent'], ascending = False))
    print("Printed 50-10")
    print("All Keys for {}: {}".format(group_in_test, len(df_all_others_together)))
    if group_in_test != 'ALPHA_BELTA':
        allowable_features_df = df_all_others_together[df_all_others_together['percent_ratio1']> \
                                                       allowable_feature_ratio]
    else:
        allowable_features_df= df_all_others_together[(df_all_others_together['group_percent']> 50 )\
                                 & (df_all_others_together['all_others_percent']< 10 )].\
                          sort_values(['group_percent'], ascending = False)
    print("Allowable keys for {} : {}".format(group_in_test, len(allowable_features_df)))
    allowable_features_df.to_csv(os.path.join(args.path, args.sample_name, args.bin_folder, 
    	'allowable_features_{}'.format(group_in_test)))


def generate_box_plot(all_others_together, all_files_in_others, group_in_test):
    group_percents = []
    all_others_percents = []
    for percent, v in all_others_together.items():
        group_percents.append(percent)
        all_others_percents.append(list(100*np.array(list(v.values()))/all_files_in_others))
    print(len(all_others_percents)   ) 
    #df_all_others_together['all_others_percent'] = df_all_others_together[
    #										'others_files_count'].div(all_files_in_others)
    #print(df_all_others_together[df_all_others_together['group_percent'] == 68.0])
    #df_all_others_together.boxplot(by ='group_percent', column =['all_others_percent'], grid = False) 
    fig = plt.figure(1)
    # Create an axes instance
    ax = fig.add_subplot(111)
    ax.boxplot(all_others_percents)
    ax.set_xticklabels(group_percents, rotation = 90)
    ax.set_xlabel("Percentage of Proteins in only {} group(%)".\
                  format(groups_display_name[group_in_test]), fontsize = 15)
    ax.set_ylabel("Respective Keys Percentage in all other groups(%)", fontsize = 15)
    ax.set_title("Keys Percentage Distribution in {} Vs all other groups(%)".\
                 format(groups_display_name[group_in_test]), fontsize = 15)
    fig = plt.gcf()
    fig.set_size_inches(18.5, 10.5)
    plt.savefig(os.path.join(args.path, args.sample_name, args.bin_folder, 'group_analysis', 
    	"keys_distribution_within_across_groups_{}.png".format(group_in_test)), bbox_inches='tight')
    plt.show()

def generate_unique_keys_per_group(df, group_in_test, groups, top_n, x_axis_step, y_axis_step):
    inter_keys = []
    for group in groups:
        if group != group_in_test:
            print(group)
            df_g = pd.read_csv(os.path.join(args.path, args.sample_name, args.bin_folder, 
            	'keys_percents_{}_0.csv'.format(group)), header = 0)
            df_inter = df[df['key'].isin(df_g['key'].values)]
            print('Intersection keys count so far:',len(inter_keys))
            inter_keys = inter_keys + list(df_inter['key'].values)
    plt.show()
    plt.close()
    print('Intersection keys in group: ', len(set(inter_keys)))
    df_unique_group = df[~df['key'].isin(list(set(inter_keys)))]
    print('Unique keys in group: ',len(df_unique_group))
    df_unique_group.to_csv(os.path.join(args.path, args.sample_name, args.bin_folder, 
    	'unique_keys_{}.csv'.format(group_in_test)))
    top_n_percents = plot_unique_keys(df_unique_group, 
                            groups_display_name[group_in_test], 
                            len(df), len(df_unique_group), top_n, x_axis_step, y_axis_step)
    
    plt.show()

    df_top_n = df_unique_group[df_unique_group['rounded_key_percent'].isin(top_n_percents)]
    # print(df_top_n)

def generate_dataset_allawable_features(groups):
    all_selected_features = []
    for group in groups:
        df_allowable_features = pd.read_csv(os.path.join(args.path, args.sample_name, args.bin_folder,
         'allowable_features_{}'.format(group)), header = 0)
        df_unique = pd.read_csv(os.path.join(args.path, args.sample_name, args.bin_folder, 
        	'unique_keys_{}.csv'.format(group)), header = 0)
        print('unique_keys_{}.csv'.format(group_in_test))
        features = df_allowable_features['key'].values
        unique = df_unique['key'].values
        print("Group: {}, Features: {}, Unique: {}".format(group, len(features), len(unique)))
        all_selected_features = list(set(all_selected_features).union(set(features)).union(set(unique)))
        print("Features so far:{}".format(len(all_selected_features)))
    print("Total dataset selected features: {}".format(len(all_selected_features)) )
    with open(os.path.join(args.path, args.sample_name, args.bin_folder, 
    	'localFeatureSelection_theta29_dist35_custom.txt'), 'w') as f:
        for item in all_selected_features:
            f.write("%s\n" % item)
    print("Selected features are saved at {}".format(os.path.join(args.path, 
    	args.sample_name, args.bin_folder, 'localFeatureSelection_theta29_dist35_custom.txt')))
    


if __name__ == '__main__':

	args = parser.parse_args()

	sample_name = args.sample_name
	if not os.path.exists(os.path.join(args.path, args.sample_name, args.bin_folder, 'group_analysis')):
		os.makedirs(os.path.join(args.path, args.sample_name, args.bin_folder, 'group_analysis'))
	groups_display_name = { "ALPHA_and_BELTA": "ALPHA+BELTA",
	                      "ALPHA_or_BELTA": "ALPHA/BELTA",
	                      "ALPHA": "ALPHA",
	                      "BELTA": "BELTA"}
	top_n = 1
	x_axis_step = {'BELTA': 0.25, 'ALPHA': 0.85, 'ALPHA_or_BELTA': 1, 'ALPHA_and_BELTA': 1}
	y_axis_step = {'BELTA': 300, 'ALPHA': 20, 'ALPHA_or_BELTA': 5000, 'ALPHA_and_BELTA': 100}


	groups = [re.search('keys_percents_(.*)_0', ntpath.basename(file)).group(1) for file in glob.glob(os.path.join(args.path, args.sample_name, args.bin_folder, 'keys_percents_*_0.csv'))]
	print(groups)
	for group in groups:
		print("#### {} #####".format(group))
		group_in_test = group		
		df = pd.read_csv(os.path.join(args.path, args.sample_name, args.bin_folder, 
			'keys_percents_{}_0.csv'.format(group_in_test)), header = 0)
		df['rounded_key_percent'] = df['key_percent'].round()

		print(df.head(5))
		print(len(df))

		#### Plot 1: Plot two-sided plot with percent distribution of keys in specific group
		all_others_together, all_files_in_others = generate_two_sided_plot(group_in_test, groups)

		#### Plot 2: Histogram type of plot?? #####
		df_all_others_together = generate_all_others_together_df(all_others_together, all_files_in_others, 8)

		generate_allowable_features(df_all_others_together, group_in_test, 1.5)

		#### Plot 3: Very Important Plot: Could be used in publication
		generate_box_plot(all_others_together, all_files_in_others, group_in_test)

		#### Plot 4: Distribution of Unique keys in specific group w.r.t to its Proteins Percentage
		generate_unique_keys_per_group(df, group_in_test, groups, top_n,                                
	                               x_axis_step[group_in_test], y_axis_step[group_in_test])

	generate_dataset_allawable_features(groups)



