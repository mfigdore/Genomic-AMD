data_fn = '/Users/mattf/Documents/Spr19/448/final/EyeGEx_retina_combined_genelevel_expectedcounts_byrid_nooutliers.counts.matrix.tsv'
data = pd.read_table(data_fn, sep = '\t', header = 0, index_col = 0)
data # genes are rows. samples are columns

cells = list(data)
print(cells, len(cells))
contained = data > 6
data = data[contained.sum(axis=1) >= 6]
data = np.transpose(data) #????
data = pd.DataFrame(scale(data))
data

pca = PCA(n_components=20) # reduce the dimensionality of the GENES
pca_samples = pca.fit_transform(data) # 10 min...
print(pca.explained_variance_ratio_)
print('SUM: ')
print(np.sum(pca.explained_variance_ratio_))
print(pca_samples.shape)