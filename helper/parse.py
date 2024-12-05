def parse_edge_list(input_file, output_file):
    with open(input_file) as f:
        n, m = map(int, f.readline().split())
        matrix = [[0 for _ in range(n)] for _ in range(n)]
        for _ in range(m):
            v1, v2 = map(int, f.readline().split())
            matrix[v1][v2] = matrix[v2][v1] = 1
            
    with open(output_file, 'w') as f:
        f.write("1\n")
        f.write(f"{n}\n")
        for row in matrix:
            f.write(' '.join(map(str, row[:n])) + '\n')

if __name__ == "__main__":
    parse_edge_list("ind0.in", "output.txt")