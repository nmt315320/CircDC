from pygraph.classes.digraph import digraph
import itertools

class MapReduce:
    __doc__ = '''提供map_reduce功能'''

    @staticmethod
    def map_reduce(i, mapper, reducer):
        """
        map_reduce方法
        :param i: 需要MapReduce的集合
        :param mapper: 自定义mapper方法
        :param reducer: 自定义reducer方法
        :return: 以自定义reducer方法的返回值为元素的一个列表
        """
        intermediate = []  # 存放所有的(intermediate_key, intermediate_value)
        for (key, value) in i.items():
            intermediate.extend(mapper(key, value))

        # sorted返回一个排序好的list，因为list中的元素是一个个的tuple，key设定按照tuple中第几个元素排序
        # groupby把迭代器中相邻的重复元素挑出来放在一起,key设定按照tuple中第几个元素为关键字来挑选重复元素
        # 下面的循环中groupby返回的key是intermediate_key，而group是个list，是1个或多个
        # 有着相同intermediate_key的(intermediate_key, intermediate_value)
        groups = {}
        for key, group in itertools.groupby(sorted(intermediate, key=lambda im: im[0]), key=lambda x: x[0]):
            groups[key] = [y for x, y in group]
        # groups是一个字典，其key为上面说到的intermediate_key，value为所有对应intermediate_key的intermediate_value
        # 组成的一个列表
        return [reducer(intermediate_key, groups[intermediate_key]) for intermediate_key in groups]


from pygraph.classes.digraph import digraph


class PRIterator:
    __doc__ = '''计算一张图中的PR值'''

    def __init__(self, dg):
        self.damping_factor = 0.85  # 阻尼系数,即α
        self.max_iterations = 100  # 最大迭代次数
        self.min_delta = 0.00001  # 确定迭代是否结束的参数,即ϵ
        self.graph = dg

    def page_rank(self):
        #  先将图中没有出链的节点改为对所有节点都有出链
        for node in self.graph.nodes():
            if len(self.graph.neighbors(node)) == 0:
                for node2 in self.graph.nodes():
                    digraph.add_edge(self.graph, (node, node2))

        nodes = self.graph.nodes()
        graph_size = len(nodes)

        if graph_size == 0:
            return {}
        page_rank = dict.fromkeys(nodes, 1.0 / graph_size)  # 给每个节点赋予初始的PR值
        damping_value = (1.0 - self.damping_factor) / graph_size  # 公式中的(1−α)/N部分

        flag = False
        for i in range(self.max_iterations):
            change = 0
            for node in nodes:
                rank = 0
                for incident_page in self.graph.incidents(node):  # 遍历所有“入射”的页面
                    rank += self.damping_factor * (page_rank[incident_page] / len(self.graph.neighbors(incident_page)))
                rank += damping_value
                change += abs(page_rank[node] - rank)  # 绝对值
                page_rank[node] = rank

            print("This is NO.%s iteration" % (i + 1))
            print(page_rank)

            if change < self.min_delta:
                flag = True
                break
        if flag:
            print("finished in %s iterations!" % node)
        else:
            print("finished out of 100 iterations!")
        return page_rank


if __name__ == '__main__':
    dg = digraph()

    dg.add_nodes(["A", "B", "C", "D", "E"])

    dg.add_edge(("A", "B"))
    dg.add_edge(("A", "C"))
    dg.add_edge(("A", "D"))
    dg.add_edge(("B", "D"))
    dg.add_edge(("C", "E"))
    dg.add_edge(("D", "E"))
    dg.add_edge(("B", "E"))
    dg.add_edge(("E", "A"))

    pr = PRIterator(dg)
    page_ranks = pr.page_rank()

    print("The final page rank is\n", page_ranks)