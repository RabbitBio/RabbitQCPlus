#include <cstdio>

struct th_ass{
    int pugz_t;
    int w_t;
    int pigz_t;
};

th_ass thread_assignment_se3[33]={
    {0, 0, 0}, // 0 0
    {0, 1, 0}, // 75 1
    {0, 1, 0}, // 75 2
    {0, 1, 2}, // 133 3
    {0, 1, 3}, // 189 4
    {2, 1, 3}, // 205 5
    {2, 1, 4}, // 260 6
    {2, 2, 4}, // 262 7
    {2, 2, 5}, // 340 8
    {2, 2, 6}, // 400 9
    {2, 2, 7}, // 460 10
    {2, 2, 7}, // 460 11
    {2, 3, 8}, // 511 12
    {2, 3, 9}, // 566 13
    {2, 3, 10}, // 621 14
    {2, 3, 11}, // 666 15
    {2, 3, 11}, // 666 16
    {3, 3, 11}, // 676 17
    {3, 3, 12}, // 680 18
    {3, 4, 12}, // 731 19
    {4, 4, 13}, // 786 20
    {4, 5, 14}, // 836 21
    {5, 5, 14}, // 841 22
    {5, 6, 15}, // 890 23
    {6, 6, 15}, // 896 24
    {6, 7, 16}, // 900 25
    {7, 7, 16}, // 948 26
    {7, 8, 17}, // 980 27
    {7, 8, 17}, // 980 28
    {8, 9, 17}, // 1000 29
    {8, 9, 18}, // 1002 30
    {8, 10, 18}, // 1002 31
    {8, 10, 18} // 1002 32
};
th_ass thread_assignment_se2[33]={
    {0, 0, 0}, // 0 0
    {0, 1, 0}, // 75 1
    {0, 1, 0}, // 75 2
    {0, 1, 2}, // 133 3
    {0, 1, 3}, // 205 4
    {0, 1, 4}, // 260 5
    {0, 2, 4}, // 262 6
    {0, 2, 5}, // 340 7
    {0, 2, 6}, // 400 8
    {0, 2, 7}, // 460 9
    {0, 2, 7}, // 460 10
    {0, 3, 8}, // 511 11
    {0, 3, 9}, // 566 12
    {0, 3, 10}, // 621 13
    {0, 3, 11}, // 676 14
    {0, 3, 12}, // 680 15
    {0, 4, 12}, // 731 16
    {0, 4, 13}, // 786 17
    {0, 4, 14}, // 841 18
    {0, 5, 15}, // 896 19
    {0, 5, 16}, // 900 20
    {0, 5, 16}, // 948 21
    {0, 6, 17}, // 980 22
    {0, 6, 17}, // 1000 23
    {0, 7, 18}, // 1050 24
    {0, 7, 19}, // 1100 25
    {0, 8, 20}, // 1150 26
    {0, 8, 20}, // 1150 27
    {0, 8, 21}, // 1150 28
    {0, 9, 22}, // 1150 29
    {0, 9, 23}, // 1150 30
    {0, 10, 23}, // 1150 31
    {0, 10, 24} // 1150 32
};
th_ass thread_assignment_se1[33]={
    {0, 0, 0}, // 0 0
    {0, 1, 0}, // 189 1
    {2, 1, 0}, // 260 2
    {2, 2, 0}, // 460 3
    {2, 3, 0}, // 666 4
    {2, 3, 0}, // 666 5
    {3, 3, 0}, // 680 6
    {4, 4, 0}, // 836 7
    {5, 4, 0}, // 890 8
    {6, 4, 0}, // 900 9
    {7, 5, 0}, // 980 10
    {7, 5, 0}, // 980 11
    {8, 6, 0}, // 1002 12
    {8, 7, 0}, // 1002 13
    {8, 8, 0}, // 1002 14
    {8, 9, 0}, // 1002 15
    {8, 10, 0}, // 1002 16
    {8, 11, 0}, // 1002 17
    {8, 12, 0}, // 1002 18
    {8, 13, 0}, // 1002 19
    {8, 14, 0}, // 1002 20
    {8, 15, 0}, // 1002 21
    {8, 16, 0}, // 1002 22
    {8, 17, 0}, // 1002 23
    {8, 18, 0}, // 1002 24
    {8, 19, 0}, // 1002 25
    {8, 20, 0}, // 1002 26
    {8, 21, 0}, // 1002 27
    {8, 22, 0}, // 1002 28
    {8, 23, 0}, // 1002 29
    {8, 24, 0}, // 1002 30
    {8, 25, 0}, // 1002 31
    {8, 26, 0} // 1002 32
};
th_ass thread_assignment_se0[33]={
    {0, 0, 0}, // 0 0
    {0, 1, 0}, // 260 1
    {0, 2, 0}, // 460 2
    {0, 3, 0}, // 680 3
    {0, 4, 0}, // 900 4
    {0, 5, 0}, // 980 5
    {0, 6, 0}, // 1150 6
    {0, 7, 0}, // 1300 7
    {0, 8, 0}, // 1300 8
    {0, 9, 0}, // 1300 9
    {0, 10, 0}, // 1300 10
    {0, 11, 0}, // 1300 11
    {0, 12, 0}, // 1300 12
    {0, 13, 0}, // 1300 13
    {0, 14, 0}, // 1300 14
    {0, 15, 0}, // 1300 15
    {0, 16, 0}, // 1300 16
    {0, 17, 0}, // 1300 17
    {0, 18, 0}, // 1300 18
    {0, 19, 0}, // 1300 19
    {0, 20, 0}, // 1300 20
    {0, 21, 0}, // 1300 21
    {0, 22, 0}, // 1300 22
    {0, 23, 0}, // 1300 23
    {0, 24, 0}, // 1300 24
    {0, 25, 0}, // 1300 25
    {0, 26, 0}, // 1300 26
    {0, 27, 0}, // 1300 27
    {0, 28, 0}, // 1300 28
    {0, 29, 0}, // 1300 29
    {0, 30, 0}, // 1300 30
    {0, 31, 0}, // 1300 31
    {0, 32, 0} // 1300 32
};
th_ass thread_assignment_pe3[65]={
    {0, 0, 0}, // 0 0
    {0, 1, 0}, // 48 1
    {0, 2, 0}, // 48 2
    {0, 1, 1}, // 90 3
    {0, 2, 1}, // 90 4
    {0, 1, 2}, // 100 5
    {2, 1, 2}, // 159 6
    {2, 2, 2}, // 159 7
    {2, 1, 3}, // 185 8
    {2, 2, 3}, // 246 9
    {2, 3, 3}, // 246 10
    {2, 2, 4}, // 314 11
    {2, 3, 4}, // 314 12
    {2, 2, 5}, // 344 13
    {2, 3, 5}, // 408 14
    {2, 4, 5}, // 408 15
    {2, 3, 6}, // 480 16
    {2, 4, 6}, // 480 17
    {2, 3, 7}, // 486 18
    {2, 4, 7}, // 552 19
    {2, 5, 7}, // 552 20
    {2, 4, 8}, // 608 21
    {2, 5, 8}, // 613 22
    {2, 6, 8}, // 613 23
    {2, 6, 9}, // 666 24
    {2, 6, 9}, // 666 25
    {3, 6, 9}, // 679 26
    {3, 6, 9}, // 679 27
    {3, 7, 10}, // 745 28
    {3, 7, 10}, // 745 29
    {3, 7, 11}, // 746 30
    {3, 7, 11}, // 750 31
    {4, 7, 11}, // 811 32
    {4, 7, 11}, // 811 33
    {4, 8, 12}, // 836 34
    {5, 8, 12}, // 877 35
    {5, 8, 12}, // 877 36
    {5, 9, 13}, // 883 37
    {5, 9, 13}, // 890 38
    {6, 9, 13}, // 940 39
    {7, 9, 13}, // 943 40
    {7, 9, 13}, // 943 41
    {7, 10, 14}, // 980 42
    {8, 10, 14}, // 1002 43
    {8, 10, 14}, // 1002 44
    {8, 10, 14}, // 1002 45
    {8, 11, 14}, // 1002 46
    {8, 11, 14}, // 1002 47
    {8, 11, 15}, // 1002 48
    {8, 12, 15}, // 1002 49
    {8, 12, 16}, // 1002 50
    {8, 12, 16}, // 1002 51
    {8, 13, 17}, // 1002 52
    {8, 13, 17}, // 1002 53
    {8, 13, 18}, // 1002 54
    {8, 14, 18}, // 1002 55
    {8, 14, 19}, // 1002 56
    {8, 14, 19}, // 1002 57
    {8, 15, 20}, // 1002 58
    {8, 15, 20}, // 1002 59
    {8, 15, 21}, // 1002 60
    {8, 16, 22}, // 1002 61
    {8, 16, 23}, // 1002 62
    {8, 16, 24}, // 1002 63
    {8, 16, 25} // 1002 64
};
th_ass thread_assignment_pe2[65]={
    {0, 0, 0}, // 0 0
    {0, 1, 0}, // 48 1
    {0, 2, 0}, // 48 2
    {0, 1, 1}, // 90 3
    {0, 2, 1}, // 90 4
    {0, 1, 2}, // 159 5
    {0, 2, 2}, // 159 6
    {0, 1, 3}, // 185 7
    {0, 2, 3}, // 246 8
    {0, 3, 3}, // 246 9
    {0, 2, 4}, // 314 10
    {0, 3, 4}, // 314 11
    {0, 2, 5}, // 344 12
    {0, 3, 5}, // 408 13
    {0, 4, 5}, // 408 14
    {0, 3, 6}, // 480 15
    {0, 4, 6}, // 480 16
    {0, 3, 7}, // 486 17
    {0, 4, 7}, // 552 18
    {0, 5, 7}, // 552 19
    {0, 4, 8}, // 608 20
    {0, 5, 8}, // 613 21
    {0, 6, 8}, // 613 22
    {0, 6, 9}, // 679 23
    {0, 6, 9}, // 679 24
    {0, 6, 10}, // 745 25
    {0, 7, 10}, // 745 26
    {0, 7, 11}, // 746 27
    {0, 7, 11}, // 811 28
    {0, 8, 11}, // 811 29
    {0, 8, 12}, // 877 30
    {0, 8, 12}, // 877 31
    {0, 9, 13}, // 883 32
    {0, 9, 13}, // 943 33
    {0, 9, 13}, // 943 34
    {0, 10, 14}, // 1009 35
    {0, 10, 14}, // 1009 36
    {0, 10, 15}, // 1021 37
    {0, 11, 15}, // 1075 38
    {0, 11, 15}, // 1075 39
    {0, 11, 16}, // 1137 40
    {0, 12, 16}, // 1137 41
    {0, 12, 17}, // 1167 42
    {0, 12, 17}, // 1200 43
    {0, 13, 17}, // 1200 44
    {0, 13, 18}, // 1260 45
    {0, 13, 18}, // 1260 46
    {0, 14, 19}, // 1266 47
    {0, 14, 19}, // 1300 48
    {0, 14, 19}, // 1300 49
    {0, 15, 20}, // 1300 50
    {0, 15, 20}, // 1300 51
    {0, 15, 20}, // 1300 52
    {0, 16, 21}, // 1300 53
    {0, 16, 22}, // 1300 54
    {0, 16, 23}, // 1300 55
    {0, 17, 24}, // 1300 56
    {0, 17, 24}, // 1300 57
    {0, 17, 25}, // 1300 58
    {0, 17, 25}, // 1300 59
    {0, 18, 25}, // 1300 60
    {0, 18, 25}, // 1300 61
    {0, 18, 26}, // 1300 62
    {0, 18, 26}, // 1300 63
    {0, 18, 26} // 1300 64
};
th_ass thread_assignment_pe1[65]={
    {0, 0, 0}, // 0 0
    {0, 1, 0}, // 100 1
    {2, 1, 0}, // 185 2
    {2, 2, 0}, // 344 3
    {2, 3, 0}, // 486 4
    {2, 4, 0}, // 608 5
    {2, 5, 0}, // 666 6
    {2, 6, 0}, // 666 7
    {3, 5, 0}, // 746 8
    {3, 6, 0}, // 750 9
    {4, 6, 0}, // 836 10
    {5, 6, 0}, // 883 11
    {5, 7, 0}, // 890 12
    {6, 7, 0}, // 940 13
    {7, 7, 0}, // 980 14
    {8, 7, 0}, // 1002 15
    {8, 8, 0}, // 1002 16
    {8, 9, 0}, // 1002 17
    {8, 10, 0}, // 1002 18
    {8, 10, 0}, // 1002 19
    {8, 11, 0}, // 1002 20
    {8, 12, 0}, // 1002 21
    {8, 13, 0}, // 1002 22
    {8, 14, 0}, // 1002 23
    {8, 15, 0}, // 1002 24
    {8, 16, 0}, // 1002 25
    {8, 17, 0}, // 1002 26
    {8, 18, 0}, // 1002 27
    {8, 19, 0}, // 1002 28
    {8, 20, 0}, // 1002 29
    {8, 21, 0}, // 1002 30
    {8, 22, 0}, // 1002 31
    {8, 23, 0}, // 1002 32
    {8, 24, 0}, // 1002 33
    {8, 25, 0}, // 1002 34
    {8, 26, 0}, // 1002 35
    {8, 27, 0}, // 1002 36
    {8, 28, 0}, // 1002 37
    {8, 29, 0}, // 1002 38
    {8, 30, 0}, // 1002 39
    {8, 31, 0}, // 1002 40
    {8, 32, 0}, // 1002 41
    {8, 33, 0}, // 1002 42
    {8, 34, 0}, // 1002 43
    {8, 35, 0}, // 1002 44
    {8, 36, 0}, // 1002 45
    {8, 37, 0}, // 1002 46
    {8, 38, 0}, // 1002 47
    {8, 39, 0}, // 1002 48
    {8, 40, 0}, // 1002 49
    {8, 41, 0}, // 1002 50
    {8, 42, 0}, // 1002 51
    {8, 43, 0}, // 1002 52
    {8, 44, 0}, // 1002 53
    {8, 45, 0}, // 1002 54
    {8, 46, 0}, // 1002 55
    {8, 47, 0}, // 1002 56
    {8, 48, 0}, // 1002 57
    {8, 49, 0}, // 1002 58
    {8, 50, 0}, // 1002 59
    {8, 51, 0}, // 1002 60
    {8, 52, 0}, // 1002 61
    {8, 53, 0}, // 1002 62
    {8, 54, 0}, // 1002 63
    {8, 55, 0} // 1002 64
};
th_ass thread_assignment_pe0[65]={
    {0, 0, 0}, // 0 0
    {0, 1, 0}, // 185 1
    {0, 2, 0}, // 344 2
    {0, 3, 0}, // 486 3
    {0, 4, 0}, // 608 4
    {0, 5, 0}, // 746 5
    {0, 6, 0}, // 883 6
    {0, 7, 0}, // 1021 7
    {0, 8, 0}, // 1167 8
    {0, 9, 0}, // 1266 9
    {0, 10, 0}, // 1300 10
    {0, 11, 0}, // 1300 11
    {0, 12, 0}, // 1300 12
    {0, 13, 0}, // 1300 13
    {0, 14, 0}, // 1300 14
    {0, 15, 0}, // 1300 15
    {0, 16, 0}, // 1300 16
    {0, 17, 0}, // 1300 17
    {0, 18, 0}, // 1300 18
    {0, 19, 0}, // 1300 19
    {0, 20, 0}, // 1300 20
    {0, 21, 0}, // 1300 21
    {0, 22, 0}, // 1300 22
    {0, 23, 0}, // 1300 23
    {0, 24, 0}, // 1300 24
    {0, 25, 0}, // 1300 25
    {0, 26, 0}, // 1300 26
    {0, 27, 0}, // 1300 27
    {0, 28, 0}, // 1300 28
    {0, 29, 0}, // 1300 29
    {0, 30, 0}, // 1300 30
    {0, 31, 0}, // 1300 31
    {0, 32, 0}, // 1300 32
    {0, 33, 0}, // 1300 33
    {0, 34, 0}, // 1300 34
    {0, 35, 0}, // 1300 35
    {0, 36, 0}, // 1300 36
    {0, 37, 0}, // 1300 37
    {0, 38, 0}, // 1300 38
    {0, 39, 0}, // 1300 39
    {0, 40, 0}, // 1300 40
    {0, 41, 0}, // 1300 41
    {0, 42, 0}, // 1300 42
    {0, 43, 0}, // 1300 43
    {0, 44, 0}, // 1300 44
    {0, 45, 0}, // 1300 45
    {0, 46, 0}, // 1300 46
    {0, 47, 0}, // 1300 47
    {0, 48, 0}, // 1300 48
    {0, 49, 0}, // 1300 49
    {0, 50, 0}, // 1300 50
    {0, 51, 0}, // 1300 51
    {0, 52, 0}, // 1300 52
    {0, 53, 0}, // 1300 53
    {0, 54, 0}, // 1300 54
    {0, 55, 0}, // 1300 55
    {0, 56, 0}, // 1300 56
    {0, 57, 0}, // 1300 57
    {0, 58, 0}, // 1300 58
    {0, 59, 0}, // 1300 59
    {0, 60, 0}, // 1300 60
    {0, 61, 0}, // 1300 61
    {0, 62, 0}, // 1300 62
    {0, 63, 0}, // 1300 63
    {0, 64, 0} // 1300 64
};



