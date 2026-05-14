from pygments.lexers.c_cpp import CppLexer
from pygments.token import Name, Keyword


class DraculaCppLexer(CppLexer):
    name = 'DraculaCpp'
    aliases = ['draculacpp']

    EXTRA_TYPES = {
        # STL types Pygments doesn't catch
        'vector', 'string', 'map', 'unordered_map', 'set', 'unordered_set',
        'pair', 'array', 'tuple', 'list', 'deque', 'optional', 'variant',
        'function', 'shared_ptr', 'unique_ptr', 'weak_ptr',
        'ifstream', 'ofstream', 'fstream', 'stringstream',
        'istringstream', 'ostringstream', 'filesystem', 'path',
        # emdw / project types
        'rcptr', 'uniqptr', 'Factor', 'DiscreteTable', 'DT', 'FProb',
        'AnyType', 'ClusterGraph', 'MessageQueue', 'RVIdType', 'RVIds',
        'RVVals', 'Idx2', 'RunConfig', 'TestCase', 'Coord',
        # common template type parameter
        'T',
    }

    def get_tokens_unprocessed(self, text):
        for index, token, value in super().get_tokens_unprocessed(text):
            if value in self.EXTRA_TYPES and token in (
                Name, Name.Other, Name.Variable, Name.Class, Name.Function
            ):
                token = Keyword.Type
            yield index, token, value
