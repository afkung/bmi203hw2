from .utils import Atom, Residue, ActiveSite

def findCentroid(points):
    """
    Compute the centroid for the vectors of a group of Active Site instance
    Input: n ActiveSite instances
    Output: the centroid vector
    """
    centroid = [0.0,0.0,0.0]
    for item in points:
        centroid = [centroid[0]+item.vector[0],centroid[1]+item.vector[1],centroid[2]+item.vector[2]]
    centroid = [centroid[0]/len(points),centroid[1]/len(points),centroid[2]/len(points)]
    return centroid

def findDistance(cluster1, cluster2):
    """
    Compute the similarity between two given groups of ActiveSite instances
    Input: two groups of ActiveSite instances
    Output: the Euclidean distane between their centroids (a floating point number)
    """
    similarity = 0.0
    
    # finding centroids
    centroid1 = findCentroid(cluster1)
    centroid2 = findCentroid(cluster2)

    #calculating Euclidean distance
    similarity = ((centroid1[0]-centroid2[0])**2 + (centroid1[1]-centroid2[1])**2 + (centroid1[2]-centroid2[2])**2)**0.5
    return similarity

def kmeansCluster(points, centroids, k):
    """
    Working function for partitioning algorithm
    Input: ActiveSite instances to cluster, centroids of cluster's vectors, number of clusters k
    Output: clustered ActiveSite instances
    """
    kmeansclusters = [[]]*k
    
    # iterating through points, finding closest centroid and assigning point to that cluster
    for item in points:
        closest_cluster = 0
        closest_distance = 99999.99
        for index in range(k):
            distance = ((item.vector[0]-centroids[index][0])**2 + (item.vector[1]-centroids[index][1])**2 + (item.vector[2]-centroids[index][2])**2)**0.5
            if distance < closest_distance:
                closest_cluster = index
                closest_distance = distance
        kmeansclusters[closest_cluster] += [item]
    return kmeansclusters

def cluster_by_partitioning(active_sites,k):
    """
    Cluster a given set of ActiveSite instances using a k-means partitioning method.

    Input: a list of ActiveSite instances, number of clusters k
    Output: a clustering of ActiveSite instances
            (this is really a list of clusters, each of which is list of
            ActiveSite instances)
    """
    centroids = [[0.0,0.0,0.0]]*k
    clusters = [[]]*k
    # initialize centroids to k random vectors    
    centroid_indices = np.random.choice(len(active_sites),k,replace = False)
    for index in range(k):
        centroids[index] = active_sites[centroid_indices[index]].vector
    # run clustering for n iterations
    iterations = 100
    for iteration in range(iterations):
        # cluster points based on centroids
        clusters = kmeansCluster(active_sites, centroids,k)
        # find new centroids
        for index in range(k):
            centroids[index] = findCentroid(clusters[index])
        print(centroids)
    return clusters

def clusterNearest(points):
    """
    Working function for hierarchical algorithm
    Input: list of ActiveSite clusters
    Output: list of ActiveSite clusters with the two closest clusters merged
    """
    nearest_distance = 99999.9
    nearest_pt1 = 0
    nearest_pt2 = 0
    for index1 in range(len(points)):
        for index2 in range(len(points)):
            if findDistance(points[index1],points[index2]) < nearest_distance and index1 != index2:
                nearest_distance = findDistance(points[index1],points[index2])
                nearest_pt1 = index1
                nearest_pt2 = index2
    points[index1] += [points[index2]]
    points.pop(index2)
    return points

def cluster_hierarchically(active_sites):
    """
    Cluster the given set of ActiveSite instances using a hierarchical algorithm.                                                                  #
    Input: a list of ActiveSite instances
    Output: a list of clusterings
            (each clustering is a list of lists of Sequence objects)
    """
    # iterate until desired number of clusters is reached
    while len(local_clusters) > number:
        local_clusters = clusterNearest(local_clusters)
        print(len(local_clusters))
    return local_clusters
    
def sum_distances(clusters):
    """
    Function that determines the tightness of clusters
    Input: list of clusters
    Output: floating point number representing average sum of distances to centroid of cluster
    """
    total_sum = 0.0
    for cluster in clusters:
        centroid = findCentroid(cluster)
        for item in cluster:
            total_sum += ((item[0]-centroid[0])**2 + (item[1]-centroid[1])**2 + (item[2]-centroid[2])**2)**0.5
    return total_sum / len(clusters)