import argparse
import logging
from collections import defaultdict

class CGCFinder:
    def __init__(self, gff_file, distance, siggenes, output, base_pair, filtered_output):
        self.gff_file = gff_file
        self.distance = distance
        self.siggenes = siggenes
        self.output = output
        self.filtered_output = filtered_output
        self.base_pair = base_pair
        self.cluster = [0, 0, 0, 0] # cazyme, tp, tf, stp
        self.num_clusters = 0
        self.contigs = defaultdict(list)
        self.load_contigs()
        self.filtered_clusters = []

    def load_contigs(self):
        with open(self.gff_file) as f:
            for line in f:
                if not line.startswith('#'):
                    row = line.rstrip().split('\t')
                    self.contigs[row[0]].append(row)
        logging.info("Contigs loaded successfully")

    def is_important(self, gene):
        return gene == 'CAZyme' or (
            gene == 'TC' and self.siggenes in ['tp', 'all', 'tp+tf', 'tp+stp']) or (
            gene == 'TF' and self.siggenes in ['tf', 'all', 'tp+tf', 'tf+stp']) or (
            gene == 'STP' and self.siggenes in ['stp', 'all', 'tp+stp', 'tf+stp'])

    def increase_cluster_count(self, gene):
        if gene == 'CAZyme':
            self.cluster[0] += 1
        elif gene == 'TC':
            self.cluster[1] += 1
        elif gene == 'TF':
            self.cluster[2] += 1
        elif gene == 'STP':
            self.cluster[3] += 1

    def is_cluster(self):
        if self.siggenes == 'all':
            return all(x > 0 for x in self.cluster)
        elif self.siggenes == 'tf':
            return self.cluster[0] > 0 and self.cluster[2] > 0
        elif self.siggenes == 'tp':
            return self.cluster[0] > 0 and self.cluster[1] > 0
        elif self.siggenes == 'stp':
            return self.cluster[0] > 0 and self.cluster[3] > 0
        elif self.siggenes == 'tp+tf':
            return self.cluster[0] > 0 and self.cluster[1] > 0 and self.cluster[2] > 0
        elif self.siggenes == 'tp+stp':
            return self.cluster[0] > 0 and self.cluster[1] > 0 and self.cluster[3] > 0
        elif self.siggenes == 'tf+stp':
            return self.cluster[0] > 0 and self.cluster[2] > 0 and self.cluster[3] > 0
        return False

    def find_near(self, contig, index):
        vals = ['null', 'null']
        k, l = index - 1, index + 1
        while k >= 0 and vals[0] == 'null':
            if self.is_important(contig[k][2]):
                vals[0] = index - k - 1
            k -= 1
        while l < len(contig) and vals[1] == 'null':
            if self.is_important(contig[l][2]):
                vals[1] = l - index - 1
            l += 1
        return vals

    def start_search(self, start_row, contig, out):
        dis, index, between, last_important = self.distance, start_row, 0, start_row
        while index < len(contig):
            index += 1
            fd = contig[index][2]
            if self.is_important(fd):
                self.increase_cluster_count(fd)
                last_important = index
                between = 0
            else:
                between += 1
            if between > dis or index >= len(contig) - 1:
                if self.is_cluster():
                    self.num_clusters += 1
                    self.write_cluster_output(contig, start_row, last_important, out)
                    self.filtered_clusters.append(contig[start_row:last_important + 1])
                self.cluster = [0, 0, 0, 0]
                return index
        return index

    def filter_clusters(self):
        filtered = []
        for cluster in self.filtered_clusters:
            valid = True
            for i in range(len(cluster) - 1):
                if int(cluster[i + 1][3]) - int(cluster[i][4]) > self.base_pair:
                    valid = False
                    break
            if valid:
                filtered.append(cluster)
        self.filtered_clusters = filtered

    def write_cluster_output(self, contig, start, end, out):
        for j in range(start, end + 1):
            fd = contig[j][2]
            if self.is_important(fd):
                up_down = self.find_near(contig, j)
                notes = contig[j][8].split(";")
                gene_id = next((note.split("=")[1] for note in notes if "ID" in note), "")
                row = [str(j), fd, str(up_down[1]), str(up_down[0]), 'CGC' + str(self.num_clusters),
                       contig[j][0], contig[j][3], contig[j][4], gene_id, contig[j][6], contig[j][8]]
            else:
                row = [str(j), 'null', 'null', 'null', 'CGC' + str(self.num_clusters), contig[j][0],
                       contig[j][3], contig[j][4], "", contig[j][6], contig[j][8]]
            out.write('\t'.join(row) + '\n')
        out.write('+++++\n')
        logging.info(f"Cluster {self.num_clusters} written to output")

    def write_filtered_output(self):
        with open(self.filtered_output, 'w') as out:
            for cluster in self.filtered_clusters:
                for j in cluster:
                    fd = j[2]
                    if self.is_important(fd):
                        up_down = self.find_near(cluster, cluster.index(j))
                        notes = j[8].split(";")
                        gene_id = next((note.split("=")[1] for note in notes if "ID" in note), "")
                        row = [str(cluster.index(j)), fd, str(up_down[1]), str(up_down[0]), 'CGC' + str(self.num_clusters),
                               j[0], j[3], j[4], gene_id, j[6], j[8]]
                    else:
                        row = [str(cluster.index(j)), 'null', 'null', 'null', 'CGC' + str(self.num_clusters), j[0],
                               j[3], j[4], "", j[6], j[8]]
                    out.write('\t'.join(row) + '\n')
                out.write('+++++\n')
        logging.info(f"Filtered clusters written to {self.filtered_output}")

    def run(self):
        with open(self.output, 'w') as out:
            for key in self.contigs:
                contig = self.contigs[key]
                self.num_clusters = 0
                i = 0
                while i < len(contig) - 1:
                    fd = contig[i][2]
                    if self.is_important(fd):
                        self.increase_cluster_count(fd)
                        i = self.start_search(i, contig, out)
                    else:
                        i += 1
        self.filter_clusters()
        self.write_filtered_output()
        logging.info(f"Processing of {self.gff_file} completed")

def main():
    parser = argparse.ArgumentParser(description='CAZyme Gene Cluster Finder')
    parser.add_argument('gff_file', help='GFF file containing genome information')
    parser.add_argument('--distance', '-d', type=int, default=2, help='The distance allowed between two signature genes')
    parser.add_argument('--siggenes', '-s', choices=['all', 'tp', 'tf', 'stp', 'tp+tf', 'tp+stp', 'tf+stp'], default='all', help='Signature genes types required.')
    parser.add_argument('--output', '-o', default='output.txt', help='Output file name')
    parser.add_argument('--filtered_output', '-f', default='filtered_output.txt', help='Filtered output file name')
    parser.add_argument('--base_pair', '-b', type=int, default=5000, help='Maximum allowed base pairs between genes in a cluster')

    args = parser.parse_args()
    
    logging.basicConfig(filename='cgc_finder.log', level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

    finder = CGCFinder(args.gff_file, args.distance, args.siggenes, args.output, args.base_pair, args.filtered_output)
    finder.run()

if __name__ == '__main__':
    main()
