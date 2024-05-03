import graphviz

def graph_dot_file(dot_file_path, output_format='png'):

    with open(dot_file_path, 'r') as f:
        dot_graph = f.read()


    graph = graphviz.Source(dot_graph)


    graph.render(filename=dot_file_path, format=output_format, cleanup=True)


    graph.view()

if __name__ == "__main__":
    dot_file_path = "greedy.dot"
    graph_dot_file(dot_file_path)
