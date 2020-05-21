import os
import sys
import math
import shutil
import pandas as pd
from scipy import stats
from collections import Counter

import matplotlib as mplt
mplt.use('Agg')
import matplotlib.pyplot as plt

import seaborn as sns
import numpy as np

class GeneUsage:
    def __init__(self, usage_dir):
        self.dfs = []
        self.individuals = []
        self.ind_map = dict()
        files = os.listdir(usage_dir)
        usage_dict = {'Memory' : 'memory', 'Background' : 'background', 'Stimulated' : 'stimulated', 'Raw' : 'raw'}
        for f in files:
            if f.find('combined_stats_usage.txt') == -1:
                continue
            ind_id = int(f.split('_')[0]) 
            self.individuals.append(ind_id)
            df = pd.read_csv(os.path.join(usage_dir, f), delim_whitespace = True)
            new_df = {'Gene' : [], 'Individual' : [], 'Usage' : [], 'UsageType' : []}
            for i in range(len(df)):
                for utype in ['Memory', 'Background', 'Stimulated', 'Raw']:
                    new_df['Gene'].append(df['Gene'][i])
                    new_df['Individual'].append(ind_id)
                    new_df['Usage'].append(df[utype][i])
                    new_df['UsageType'].append(usage_dict[utype])
            self.dfs.append(pd.DataFrame(new_df))
            self.ind_map[self.individuals[-1]] = len(self.individuals) - 1
        self.df = pd.concat(self.dfs)

    def GetIndividualUsage(self, ind_id, gene, usage_type):
        df = self.dfs[self.ind_map[ind_id]]
        return df.loc[(df['Gene'] == gene + '*01') & (df['UsageType'] == usage_type)]['Usage'].iloc[0]

    def GetDataFrameByGene(self, gene):
        return self.df.loc[self.df['Gene'] == gene + '*01'].reset_index()
            
class GenotypeMatrix:
    def __init__(self, gene, snp_txt):
        self.gene = gene
        snp_df = pd.read_csv(snp_txt, delim_whitespace = True)
        self._CheckPositions(snp_df)
        self.df = {'Position' : [], 'State' : [], 'Individual' : [], 'Gene' : [], 'Genotype' : []}
        self.genotype_dict = dict()
        for i in range(len(snp_df)):
            if snp_df['Individual'][i] in self.bad_individuals:
                continue
            snps = snp_df['Haplotype'][i].split(',')
            genotype_id = snp_df['HaplotypeId'][i].replace('H', 'G')
            if genotype_id not in self.genotype_dict:
                self.genotype_dict[genotype_id] = []
            self.genotype_dict[genotype_id].append(snp_df['Individual'][i])
            for snp in snps:
                splits = snp.split(':')
                if int(splits[0]) < 20:
                    continue
                self.df['Position'].append(int(splits[0]) + 1)
                self.df['State'].append(splits[1])
                self.df['Individual'].append(snp_df['Individual'][i])
                self.df['Gene'].append(gene + '*01')
                self.df['Genotype'].append(snp_df['Haplotype'][i])
        self.df = pd.DataFrame(self.df)
        self.positions = set(self.df['Position'])

    def _CheckPositions(self, all_snp_df):
#        print '==== ' + self.gene
        self.bad_individuals = []
        haplotype_counter = Counter(all_snp_df['HaplotypeId'])
        for i in range(len(all_snp_df)):
            perc_genotype = round(float(haplotype_counter[all_snp_df['HaplotypeId'][i]]) / len(all_snp_df) * 100)
            if round(float(haplotype_counter[all_snp_df['HaplotypeId'][i]]) / 204 * 100) < 5:
                self.bad_individuals.append(all_snp_df['Individual'][i])
#        print 'Bad individuals: ' + str(self.bad_individuals)
 
    def Gene(self):
        return self.gene

    def PosIter(self):
        for pos in self.positions:
            yield pos

    def IndIter(self):
        for ind in self.individuals:
            yield ind

    def GetDF(self):
        return self.df

    def GenotypeIter(self):
        for g in sorted(self.genotype_dict):
            yield g, self.genotype_dict[g]

    def BadIndividuals(self):
        return self.bad_individuals

def CreateGenotypeMatrixes(haplotype_dir):
    genotype_map = dict()
    files = os.listdir(haplotype_dir)
    for f in files:
        gene = f.split('.')[0]
        fname = os.path.join(haplotype_dir, f)
        genotype_map[gene] = GenotypeMatrix(gene, fname)
    return genotype_map

def OutputGenotype(genotype_map, individuals, output_dir):
    v_genes = [v for v in sorted(genotype_map) if v not in ['IGHV1-25', 'IGHV1-37', 'IGHV1-39']]
    matrix = []
    for ind in individuals:
        matrix.append([16] * len(v_genes))
    for v in v_genes:
        for g, ind_list in genotype_map[v].GenotypeIter():
            for ind in ind_list:
                matrix[individuals.index(ind)][v_genes.index(v)] = int(g[1 : ])
    sns.clustermap(matrix, xticklabels = [v[3 : ] for v in v_genes], yticklabels = [], vmax = 20, vmin = 1, cmap = 'tab20')
    plt.savefig(os.path.join(output_dir, 'genotypes.pdf'))
    plt.clf()

def OutputGenotypeMatrix(gene, gene_usage, utralong_df, genotype):
    usage_matrix = []
    usage_annot_matrix = []
    ul_matrix = []
    ul_annot_matrix = []
    genotypes = []
    for g, ind_list in genotype.GenotypeIter():
        usages = [gene_usage.GetIndividualUsage(ind_id, gene, 'raw') for ind_id in ind_list]
        usage_matrix.append([np.mean(usages)])
        usage_annot_matrix.append(["{0:.2f}".format(np.mean(usages))])
        ul_fracts = [utralong_df.loc[utralong_df['Individual'] == ind_id]['FractionUL'].iloc[0] for ind_id in ind_list]
        ul_matrix.append([np.mean(ul_fracts)])
        ul_annot_matrix.append(["{0:.2f}".format(np.mean(ul_fracts))])
        genotypes.append(g)
    fig, axes = plt.subplots(ncols = 2)
    sns.heatmap(usage_matrix, annot = np.array(usage_annot_matrix), yticklabels = genotypes, xticklabels = [],
                cmap = 'coolwarm', ax = axes[0], fmt = '', cbar = False)
    plt.sca(axes[0])
    plt.xlabel('average usage')
    sns.heatmap(ul_matrix, annot = np.array(ul_annot_matrix), yticklabels = [], xticklabels = [], cmap = 'coolwarm',
                ax = axes[1], fmt = '', cbar = False)
    plt.yticks(rotation = 0)
    plt.sca(axes[1])
    plt.xlabel('% of ultralong CDR3s')
    plt.suptitle(gene)
    plt.savefig(gene + '_usage_matrix.pdf')
    plt.clf()

def OutputSNPMatrix(gene, gene_usage, utralong_df, genotype):
    merged_df = pd.merge(gene_usage.GetDataFrameByGene(gene), genotype.GetDF(), how = 'left', on = ['Individual', 'Gene'])
    merged_df = pd.merge(merged_df, utralong_df, how = 'left', on = ['Individual'])
    for pos in sorted(genotype.PosIter()):
        pos_df = merged_df.loc[(merged_df['Position'] == pos) & (merged_df['UsageType'] == 'raw')].reset_index()
        if len(set(pos_df['State'])) == 1:
            continue
        fig, axes = plt.subplots(nrows = 2)
        state_usage_dict = dict()
        ul_frac_dict = dict()
        for i in range(len(pos_df)):
            state = pos_df['State'][i]
            if state not in state_usage_dict:
                state_usage_dict[state] = []
                ul_frac_dict[state] = []
            state_usage_dict[state].append(pos_df['Usage'][i])
            ul_frac_dict[state].append(pos_df['FractionUL'][i])
        states = sorted(state_usage_dict.keys())
        usage_matrix = [[]]
        usage_annot_matrix = [[]]
        ul_matrix = [[]]
        ul_annot_matrix = [[]]        
        for state in states:
            usage_matrix[0].append(np.mean(state_usage_dict[state]))
            usage_annot_matrix[0].append("{0:.2f}".format(np.mean(state_usage_dict[state])))
            ul_matrix[0].append(np.mean(ul_frac_dict[state]))
            ul_annot_matrix[0].append("{0:.2f}".format(np.mean(ul_frac_dict[state])))
        fig, axes = plt.subplots(nrows = 2)
        sns.heatmap(usage_matrix, annot = np.array(usage_annot_matrix), yticklabels = [], xticklabels = states,
                    cmap = 'coolwarm', ax = axes[0], fmt = '', cbar = False)
        plt.sca(axes[0])
        plt.ylabel('average usage')
        sns.heatmap(ul_matrix, annot = np.array(ul_annot_matrix), yticklabels = [], xticklabels = [], cmap = 'coolwarm',
                    ax = axes[1], fmt = '', cbar = False)
        plt.sca(axes[1])
        plt.ylabel('% of ultralong CDR3s')
        plt.suptitle(gene + ', position ' + str(pos))
        plt.savefig(gene + '_' + str(pos) + '_heatmap.pdf')
        plt.clf()

def TestUsagePerSNP(gene, gene_usage, genotype, output_dir):
    merged_df = pd.merge(gene_usage.GetDataFrameByGene(gene), genotype.GetDF(), how = 'left', on = ['Individual', 'Gene'])
    p_df = {'Gene' : [], 'Position' : [], 'P-value' : [], 'Likelihood' : [], 'UsageType' : [], 'F-value' : []}
    for pos in sorted(genotype.PosIter()):
        pos_df = merged_df.loc[(merged_df['Position'] == pos) & (merged_df['UsageType'] == 'raw')].reset_index()
        if len(set(pos_df['State'])) == 1:
            continue
        basename = gene + '_position_' + str(pos)
        states = sorted(list(set(pos_df['State'])))
        plt.figure()
        sns.boxplot(x = 'State', y = 'Usage', data = pos_df, order = states)
        plt.xlabel('')
        plt.ylabel('usage')
        plt.title(gene + ', position ' + str(pos))
        plt.savefig(os.path.join(output_dir, basename + '_usage.pdf'))
        plt.clf()
        #
        pos_df = merged_df.loc[(merged_df['Position'] == pos) & (merged_df['UsageType'] != 'raw')].reset_index()
        plt.figure()
        sns.boxplot(x = 'UsageType', y = 'Usage', hue = 'State', data = pos_df,
                    order = ['background', 'memory', 'stimulated'], hue_order = states)
        plt.title(gene + ', position ' + str(pos))
        plt.xlabel('')
        plt.ylabel('usage')
        plt.savefig(os.path.join(output_dir, basename + '_per_usage_type.pdf'))
        plt.clf()
        pos_df.to_csv(os.path.join(output_dir, basename + '_per_usage_type.csv'), sep = ',', index = False)
        #
        pos_df = merged_df.loc[(merged_df['Position'] == pos)]
        states = list(set(pos_df['State']))
        utypes = ['background', 'memory', 'stimulated', 'raw'] #usage_types:
        for utype in utypes:
            usages = [pos_df.loc[(pos_df['State'] == state) & (pos_df['UsageType'] == utype)]['Usage'] for state in states]
            pi_stats = ('NA', 1)
            if len(usages) == 2:
                pi_stats = stats.f_oneway(usages[0], usages[1])
            else:
                pi_stats = stats.f_oneway(usages[0], usages[1], usages[2])
            p_df['Gene'].append(gene[3 : ])
            p_df['Position'].append(pos)
            p_df['P-value'].append(pi_stats[1])
            p_df['F-value'].append(pi_stats[0])
            p_df['Likelihood'].append(-math.log(pi_stats[1], 10))
            p_df['UsageType'].append(utype)
    return pd.DataFrame(p_df)

def TestUsagePerGene(gene, gene_usage, output_dir):
    gene_usage_df = gene_usage.GetDataFrameByGene(gene)
    plt.figure()
    sns.boxplot(x = 'UsageType', y = 'Usage', data = gene_usage_df, order = ['background', 'memory', 'stimulated'],
                palette = 'Set1')
    plt.xlabel('')
    plt.ylabel('usage')
    plt.title(gene)
    plt.savefig(os.path.join(output_dir, gene + '_dynamics.pdf'))
    pi_stats = stats.f_oneway(gene_usage_df.loc[gene_usage_df['UsageType'] == 'background']['Usage'],
                              gene_usage_df.loc[gene_usage_df['UsageType'] == 'memory']['Usage'],
                              gene_usage_df.loc[gene_usage_df['UsageType'] == 'stimulated']['Usage'])
    return gene, pi_stats

def TestULCDR3s(gene, ul_cdr3s, genotype, output_dir):
    merged_df = pd.merge(ul_cdr3s, genotype.GetDF(), how = 'left', on = ['Individual'])
    p_df = {'Gene' : [], 'Position' : [], 'P-value' : [], 'Likelihood' : [], 'F-value' : []}
    for pos in genotype.PosIter():
        pos_df = merged_df.loc[(merged_df['Position'] == pos) & (merged_df['Gene'] == gene + '*01')]
        if len(set(pos_df['State'])) == 1:
            continue
        plt.figure()
        sns.boxplot(x = 'State', y = 'FractionUL', data = pos_df, order = sorted(list(set(pos_df['State']))))
        plt.ylabel('% ultralong CDR3s')
        plt.title(gene + ', position ' + str(pos))
        plt.savefig(os.path.join(output_dir, gene + '_position_' + str(pos) + '_ul_cdr3s.pdf'))
        plt.clf()
        states = list(set(pos_df['State']))
        usages = [pos_df.loc[pos_df['State'] == state]['FractionUL'] for state in states]
        pi_stats = (0, 1)
        if len(usages) == 2:
            pi_stats = stats.f_oneway(usages[0], usages[1])
        else:
            pi_stats = stats.f_oneway(usages[0], usages[1], usages[2])
        p_df['Gene'].append(gene[3 : ])
        p_df['Position'].append(pos)
        p_df['P-value'].append(pi_stats[1])
        p_df['F-value'].append(pi_stats[0])
        p_df['Likelihood'].append(-math.log(pi_stats[1], 10))
    return pd.DataFrame(p_df) 

def ReadULCDRTxt(txt):
    df = {'Individual' : [], 'FractionUL' : []}
    ul_df = pd.read_csv(txt, delim_whitespace = True)
    for i in range(len(ul_df)):
        df['Individual'].append(int(ul_df['Individual'][i]))
        ul_perc = [ul_df[str(perc)][i] for perc in range(150, 234, 10)]
        df['FractionUL'].append(sum(ul_perc))
    return pd.DataFrame(df)

def ReadTiterDF(titer_txt):
    titer_df = pd.read_csv(titer_txt, delim_whitespace = True)
    exp_df = {'Individual' : [], 'PreVacc' : [], 'InitialVacc' : [], 'BoosterVacc' : [], 'BoosterResponseVacc' : []}
    for i in range(len(titer_df)):
        exp_df['Individual'].append(titer_df['Individual'][i])
        exp_df['PreVacc'].append(math.pow(2, titer_df['PreVacc'][i]))
        exp_df['InitialVacc'].append(math.pow(2, titer_df['InitialVacc'][i]))
        exp_df['BoosterVacc'].append(math.pow(2, titer_df['BoosterVacc'][i]))
        exp_df['BoosterResponseVacc'].append(math.pow(2, titer_df['BoosterResponseVacc'][i]))
    return titer_df #pd.DataFrame(exp_df)

def TestAntibodyTiter(gene, titer_df, genotype, output_dir):
    merged_df = pd.merge(titer_df, genotype.GetDF(), how = 'left', on = ['Individual'])
    p_df = {'Gene' : [], 'Position' : [], 'P-value' : [], 'Likelihood' : [], 'TimePoint' : [], 'F-value' : []}
    for pos in genotype.PosIter():
        pos_df = merged_df.loc[(merged_df['Position'] == pos) & (merged_df['Gene'] == gene + '*01')]
        if len(set(pos_df['State'])) == 1:
            continue
        time_points = ['PreVacc', 'InitialVacc', 'BoosterVacc', 'BoosterResponseVacc']
        basename = gene + '_position_' + str(pos) + '_titer'
        for tp in time_points:
            sns.boxplot(x = 'State', y = tp, data = pos_df, order = sorted(list(set(pos_df['State']))))
            plt.ylabel('antibody titer')
            plt.title(gene + ', position ' + str(pos))
            plt.savefig(os.path.join(output_dir, basename + '_' + tp.lower() + '.pdf'))
            plt.clf()
            states = list(set(pos_df['State']))
            usages = [pos_df.loc[pos_df['State'] == state][tp] for state in states]
            pi_stats = ('NA', 1)
            if len(usages) == 2:
                pi_stats = stats.kruskal(usages[0], usages[1])
            else:
                pi_stats = stats.kruskal(usages[0], usages[1], usages[2])
            p_df['Gene'].append(gene[3 : ])
            p_df['Position'].append(pos)
            p_df['P-value'].append(pi_stats[1])
            p_df['F-value'].append(pi_stats[0])
            p_df['Likelihood'].append(-math.log(pi_stats[1], 10))
            p_df['TimePoint'].append(tp)
    return pd.DataFrame(p_df)

def PrepareOutputDir(output_dir):
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
    os.mkdir(output_dir)

def OutputUsageStats(usage_dfs, output_dir):
    usage_df = pd.concat(usage_dfs)
    raw_df = usage_df.loc[usage_df['UsageType'] == 'raw']
    sns.swarmplot(x = 'Gene', y = 'Likelihood', data = raw_df)
    plt.ylabel('$-log_{10}$(p-value)')
    plt.savefig(os.path.join(output_dir, 'usage_per_snp_manhattan.pdf'))
    plt.clf()
    raw_df.to_csv(os.path.join(output_dir, 'usage_per_snp.csv'), sep = ',', index = False)

def OutputUltralongStats(ul_cdr3_dfs, output_dir):
    ul_cdr3_df = pd.concat(ul_cdr3_dfs)
    plt.figure()
    sns.swarmplot(x = 'Gene', y = 'Likelihood', data = ul_cdr3_df)
    plt.ylabel('$-log_{10}$(p-value)')
    plt.savefig(os.path.join(output_dir, 'ul_cdr3s_per_snp_manhattan.pdf'))
    plt.clf()
    ul_cdr3_df.to_csv(os.path.join(output_dir, 'ul_cdr3s_per_snp.csv'), sep = ',', index = False)

def OutputTiterStats(titer_dfs, output_dir):
    titer_df = pd.concat(titer_dfs)
    plt.figure()
    sns.swarmplot(x = 'Gene', y = 'Likelihood', data = titer_df)
    plt.ylabel('$-log_{10}$(p-value)')
    plt.savefig(os.path.join(output_dir, 'titer_per_snp_manhattan.pdf'))
    plt.clf()
    titer_df.to_csv(os.path.join(output_dir, 'titer_per_snp.csv'), sep = ',', index = False)

def OutputUsageDynamicsStats(usage_dynamics_stats, output_dir):
    df = {'Gene' : [], 'P-value' : [], 'F-value' : [], 'Likelihood' : []}
    for gene, stats in usage_dynamics_stats:
        df['Gene'].append(gene[3 : ])
        df['P-value'].append(stats[0])
        df['F-value'].append(stats[1])
        df['Likelihood'].append(-math.log(stats[1], 10))
    df = pd.DataFrame(df)
    plt.figure()
    sns.barplot(x = 'Gene', y = 'Likelihood', data = df)
    plt.ylabel('$-log_{10}$(p-value)')
    plt.savefig(os.path.join(output_dir, 'usage_dynamics_manhattan.pdf'))
    plt.clf()
    df.to_csv(os.path.join(output_dir, 'usage_dynamics.csv'), sep = ',', index = False)

def main(haplotype_dir, usage_dir, ul_cdr3s_txt, antibody_titer_txt, output_dir):
    PrepareOutputDir(output_dir)
    genotype_map = CreateGenotypeMatrixes(haplotype_dir)
    gene_usage = GeneUsage(usage_dir)
    ul_cdr3_df = ReadULCDRTxt(ul_cdr3s_txt)
    titer_df = ReadTiterDF(antibody_titer_txt)
    OutputGenotype(genotype_map, sorted(list(ul_cdr3_df['Individual'])), output_dir)
    usage_dfs = []
    usage_dynamics_stats = []
    ul_cdr3_dfs = []
    titer_dfs = []
    for gene in sorted(genotype_map):
        if gene not in ['IGHV1-7', 'IGHV1-10', 'IGHV1-14', 'IGHV1-17', 'IGHV1-20', 'IGHV1-21', 'IGHV1-27', 'IGHV1-30']:
            continue
        print 'Processing ' + gene + '...'
        usage_dfs.append(TestUsagePerSNP(gene, gene_usage, genotype_map[gene], output_dir))
        usage_dynamics_stats.append(TestUsagePerGene(gene, gene_usage, output_dir))
        ul_cdr3_dfs.append(TestULCDR3s(gene, ul_cdr3_df, genotype_map[gene], output_dir))
        titer_dfs.append(TestAntibodyTiter(gene, titer_df, genotype_map[gene], output_dir))
    OutputUsageStats(usage_dfs, output_dir)
    OutputUltralongStats(ul_cdr3_dfs, output_dir)
    OutputTiterStats(titer_dfs, output_dir)
    OutputUsageDynamicsStats(usage_dynamics_stats, output_dir)

if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
