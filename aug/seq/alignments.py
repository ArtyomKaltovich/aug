import sys
import math
from abc import abstractmethod

REL_TOL = 0.0
ABS_TOL = 1e-6


def are_equal(a, b):
    if type(a) == float or type(b) == float:
        return math.isclose(a, b, rel_tol=REL_TOL, abs_tol=ABS_TOL)
    else:
        return a == b


class BaseAlignment:
    def init_distance_matrix(self, seq1, seq2):
        pass

    @abstractmethod
    def calculate_distance(self, seq1, seq2, distances, i, j):
        pass

    @abstractmethod
    def score(self, distances):
        pass

    @abstractmethod
    def reconstruct_answer(self, seq1, seq2, distance, swap_case_on_mismatch=True):
        result1, result2 = self._reconstruct_answer(seq1, seq2, distance, swap_case_on_mismatch=True)
        if swap_case_on_mismatch:
            self._swap_case_on_mismatch(result1, result2)
        return "".join(reversed(result1)), "".join(reversed(result2))

    @abstractmethod
    def _reconstruct_answer(self, seq1, seq2, distance, swap_case_on_mismatch=True):
        pass

    def _swap_case_on_mismatch(self, result1, result2):
        for i, (a, b) in enumerate(zip(result1, result2)):
            if a != b:
                result2[i] = b.swapcase()


class Levinshtein(BaseAlignment):
    def __init__(self):
        super().__init__()
        self.match_score = 0
        self.gap_score = -1
        self.mismatch_score = -1

    def init_distance_matrix(self, seq1, seq2):
        distances = [[0] * (len(seq1) + 1) for _ in range(len(seq2) + 1)]
        for x in range(1, len(distances[0])):
            distances[0][x] = distances[0][x - 1] + self.gap_score
        for x in range(1, len(distances)):
            distances[x][0] = distances[x - 1][0] + self.gap_score
        return distances

    def calculate_distance(self, seq1, seq2, distances, i, j):
        match_or_mismatch_score = self._calculate_match_mismatch_score(seq1, seq2, i, j)
        i_minus_1 = i - 1
        j_minus_1 = j - 1
        distances[i][j] = max(distances[i_minus_1][j_minus_1] + match_or_mismatch_score,
                              distances[i_minus_1][j] + self.gap_score,
                              distances[i][j_minus_1] + self.gap_score)

    def score(self, distances):
        return -distances[-1][-1]

    def _reconstruct_answer(self, seq1, seq2, distance, swap_case_on_mismatch=True):
        result1 = []
        result2 = []
        self._reconstruct_until_zero_row_or_column(distance, result1, result2, seq1, seq2, swap_case_on_mismatch)
        return result1, result2

    def _calculate_match_mismatch_score(self, seq1, seq2, i, j):
        match_or_mismatch_score = self.mismatch_score if seq2[i - 1] != seq1[j - 1] else self.match_score
        return match_or_mismatch_score

    def _reconstruct_until_zero_row_or_column(self, distance, result1, result2, seq1, seq2, swap_case_on_mismatch):
        i = len(distance) - 1
        j = len(distance[0]) - 1
        while i or j:
            i_minus_1 = i - 1
            j_minus_1 = j - 1
            if self._mis_or_match_case(seq1, seq2, distance, i, j, i_minus_1, j_minus_1):
                i, j = self._on_mis_or_match(i, i_minus_1, j, j_minus_1, result1, result2, seq1, seq2)
            elif self._delete_case(seq1, seq2, distance, i, j, i_minus_1, j_minus_1):
                j = self._on_delete(j, j_minus_1, result1, result2, seq1)
            elif self._insert_case(seq1, seq2, distance, i, j, i_minus_1, j_minus_1):
                i = self._on_insert(i, i_minus_1, result1, result2, seq2)

    def _on_insert(self, i, i_minus_1, result1, result2, seq2):
        result1.append("-")
        result2.append(seq2[i_minus_1])
        i = i_minus_1
        return i

    def _on_delete(self, j, j_minus_1, result1, result2, seq1):
        result1.append(seq1[j_minus_1])
        result2.append("-")
        j = j_minus_1
        return j

    def _on_mis_or_match(self, i, i_minus_1, j, j_minus_1, result1, result2, seq1, seq2):
        result1.append(seq1[j_minus_1])
        result2.append(seq2[i_minus_1])
        i = i_minus_1
        j = j_minus_1
        return i, j

    def _mis_or_match_case(self, seq1, seq2, distance, i, j, i_minus_1, j_minus_1):
        return i and j and are_equal(distance[i_minus_1][j_minus_1],
                                     distance[i][j] - self._calculate_match_mismatch_score(seq1, seq2, i, j))

    def _insert_case(self, seq1, seq2, distance, i, j, i_minus_1, j_minus_1):
        return i and are_equal(distance[i_minus_1][j], distance[i][j] - self.gap_score)

    def _delete_case(self, seq1, seq2, distance, i, j, i_minus_1, j_minus_1):
        return j and are_equal(distance[i][j_minus_1], distance[i][j] - self.gap_score)


class NeedlemanWunsch(Levinshtein):
    def __init__(self, match_score=1, mismatch_score=-1, gap_score=-1, score_matrix=None, gap_start=None):
        self.match_score = match_score
        self.mismatch_score = mismatch_score
        self.gap_score = gap_score
        self.score_matrix = score_matrix
        # the following field are for affine gap
        self.gap_start = gap_start
        if self.gap_start is not None:
            self.min_ = -sys.maxsize - 1 if all(type(x) == int for x in
                                                (gap_start, gap_score, match_score, mismatch_score)) else float("-inf")
        self._gap_matrix_m = None
        self._gap_matrix_x = None
        self._gap_matrix_y = None

    def init_distance_matrix(self, seq1, seq2):
        if self.gap_start is not None:
            self._gap_matrix_m = [[self._init_m(i, j) for j in range(len(seq1) + 1)] for i in range(len(seq2) + 1)]
            self._gap_matrix_x = [[self._init_x(i, j) for j in range(len(seq1) + 1)] for i in range(len(seq2) + 1)]
            self._gap_matrix_y = [[self._init_y(i, j) for j in range(len(seq1) + 1)] for i in range(len(seq2) + 1)]
        else:
            return super().init_distance_matrix(seq1, seq2)

    def calculate_distance(self, seq1, seq2, distances, i, j):
        match_or_mismatch_score = self._calculate_match_mismatch_score(seq1, seq2, i, j)
        i_minus_1 = i - 1
        j_minus_1 = j - 1
        if self.gap_start is not None:
            self._calculate_distance_with_affine_gap(match_or_mismatch_score, i, j, i_minus_1, j_minus_1)
        else:
            super().calculate_distance(seq1, seq2, distances, i, j)

    def _init_m(self, i, j):
        return self.min_ if i > 0 and j == 0 or j > 0 and i == 0 else 0

    def _init_x(self, i, j):
        if i > 0 and j == 0:
            return self.min_
        elif i == 0 and j > 0:
            return self.gap_start + self.gap_score * j
        return 0

    def _init_y(self, i, j):
        if i > 0 and j == 0:
            return self.gap_start + self.gap_score * i
        elif i == 0 and j > 0:
            return self.min_
        return 0

    def _calculate_distance_with_affine_gap(self, match_or_mismatch_score, i, j, i_minus_1, j_minus_1):
        self._gap_matrix_m[i][j] = match_or_mismatch_score + max(self._gap_matrix_m[i_minus_1][j_minus_1],
                                                                 self._gap_matrix_x[i_minus_1][j_minus_1],
                                                                 self._gap_matrix_y[i_minus_1][j_minus_1])
        gap_start_score = self.gap_start + self.gap_score
        self._gap_matrix_x[i][j] = max((gap_start_score + self._gap_matrix_m[i][j_minus_1]),
                                       (self.gap_score + self._gap_matrix_x[i][j_minus_1]),
                                       (gap_start_score + self._gap_matrix_y[i][j_minus_1]))
        self._gap_matrix_y[i][j] = max((gap_start_score + self._gap_matrix_m[i_minus_1][j]),
                                       (gap_start_score + self._gap_matrix_x[i_minus_1][j]),
                                       (self.gap_score + self._gap_matrix_y[i_minus_1][j]))

    def score(self, distances):
        if self.gap_start is not None:
            return max(self._gap_matrix_m[-1][-1], self._gap_matrix_x[-1][-1], self._gap_matrix_y[-1][-1])
        else:
            return -super().score(distances)

    def _reconstruct_answer(self, seq1, seq2, distance, swap_case_on_mismatch=True):
        if self.gap_start is not None:
            result1 = []
            result2 = []
            i = len(seq2)
            j = len(seq1)
            prev_i = i
            prev_j = j
            gap_start_score = self.gap_start + self.gap_score
            score = self.score(distance)
            if self._gap_matrix_m[-1][-1] == score:
                current_matrix = "m"
            elif self._gap_matrix_x[-1][-1] == score:
                current_matrix = "x"
            elif self._gap_matrix_y[-1][-1] == score:
                current_matrix = "y"
            score_m = None
            score_x = None
            score_y = None
            while i or j:
                i_minus_1 = i - 1
                j_minus_1 = j - 1
                #print("".join(reversed(result1)), "".join(reversed(result2)))
                if self._affine_mis_or_match_case(current_matrix, i, j, prev_i, prev_j, score_m):
                    match_or_mismatch_score = self._calculate_match_mismatch_score(seq1, seq2, i, j)
                    prev_i = i
                    prev_j = j
                    i, j = super()._on_mis_or_match(i, i_minus_1, j, j_minus_1, result1, result2, seq1, seq2)
                    current_matrix = self._gap_matrix_m
                    score_m = match_or_mismatch_score
                    score_x = match_or_mismatch_score
                    score_y = match_or_mismatch_score
                elif self._affine_delete_case(current_matrix, i, j, prev_i, prev_j, score_x):
                    prev_i = i
                    prev_j = j
                    j = super()._on_delete(j, j_minus_1, result1, result2, seq1)
                    current_matrix = self._gap_matrix_x
                    score_m = gap_start_score
                    score_x = self.gap_score
                    score_y = gap_start_score
                elif self._affine_insert_case(current_matrix, i, j, prev_i, prev_j, score_y):
                    prev_i = i
                    prev_j = j
                    i = super()._on_insert(i, i_minus_1, result1, result2, seq2)
                    current_matrix = self._gap_matrix_y
                    score_m = gap_start_score
                    score_x = gap_start_score
                    score_y = self.gap_score
            return result1, result2
        else:
            return super()._reconstruct_answer(seq1, seq2, distance, swap_case_on_mismatch)

    def _affine_mis_or_match_case(self, current_matrix, i, j, prev_i, prev_j, score_m):
        return current_matrix == "m" or \
               type(current_matrix) != str and \
               are_equal(current_matrix[prev_i][prev_j], self._gap_matrix_m[i][j] + score_m)

    def _affine_delete_case(self, current_matrix, i, j, prev_i, prev_j, score_x):
        return current_matrix == "x" or \
               type(current_matrix) != str and \
               are_equal(current_matrix[prev_i][prev_j], self._gap_matrix_x[i][j] + score_x)

    def _affine_insert_case(self, current_matrix, i, j, prev_i, prev_j, score_y):
        return current_matrix == "y" or \
               type(current_matrix) != str and \
               are_equal(current_matrix[prev_i][prev_j], self._gap_matrix_y[i][j] + score_y)
