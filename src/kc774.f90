            this%wimp_lin(1) = 0.1235_dp

            this%butcher_imp(:,:) = reshape([0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,        &
            &    0.0_dp, 0.0_dp, 0.0_dp, 0.1235_dp, 0.1235_dp, 0.0_dp, 0.0_dp,      &
            &    0.0_dp, 0.0_dp, 0.0_dp, 624185399699.0_dp/4186980696204.0_dp,      &
            &    624185399699.0_dp/4186980696204.0_dp, 0.1235_dp, 0.0_dp, 0.0_dp,   &
            &    0.0_dp, 0.0_dp, 1258591069120.0_dp/10082082980243.0_dp,            &
            &    1258591069120.0_dp/10082082980243.0_dp, -322722984531.0_dp/        &
            &    8455138723562.0_dp, 0.1235_dp, 0.0_dp, 0.0_dp, 0.0_dp,             &
            &    -436103496990.0_dp/5971407786587.0_dp, -436103496990.0_dp/         &
            &    5971407786587.0_dp, -2689175662187.0_dp/11046760208243.0_dp,       &
            &    4431412449334.0_dp/12995360898505.0_dp, 0.1235_dp, 0.0_dp, 0.0_dp, &
            &    -2207373168298.0_dp/14430576638973.0_dp, -2207373168298.0_dp/      &
            &    14430576638973.0_dp, 242511121179.0_dp/3358618340039.0_dp,         &
            &    3145666661981.0_dp/7780404714551.0_dp,  5882073923981.0_dp/        &
            &    14490790706663.0_dp, 0.1235_dp, 0.0_dp, 0.0_dp, 0.0_dp,            &
            &    9164257142617.0_dp/17756377923965.0_dp, -10812980402763.0_dp/      &
            &    74029279521829.0_dp, 1335994250573.0_dp/5691609445217.0_dp,        &
            &    2273837961795.0_dp/8368240463276.0_dp, 0.1235_dp], [7,7], order=[2,1])

            this%butcher_ass_imp(:) = [0.0_dp, 0.0_dp, 9164257142617.0_dp/        &
            &                          17756377923965.0_dp, -10812980402763.0_dp/ &
            &                          74029279521829.0_dp, 1335994250573.0_dp/   &
            &                          5691609445217.0_dp,  2273837961795.0_dp/   &
            &                          8368240463276.0_dp, 0.1235_dp]

            this%butcher_exp(:,:) = reshape([ 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,      &
            &    0.0_dp, 0.0_dp, 0.0_dp, 0.247_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
            &    0.0_dp, 0.0_dp, 0.06175_dp, 2694949928731.0_dp/7487940209513.0_dp,&
            &    0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 464650059369.0_dp/        &
            &    8764239774964.0_dp, 878889893998.0_dp/2444806327765.0_dp,         &
            &    -952945855348.0_dp/12294611323341.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,   &
            &    0.0_dp, 476636172619.0_dp/8159180917465.0_dp, -1271469283451.0_dp/&
            &    7793814740893.0_dp, -859560642026.0_dp/4356155882851.0_dp,        &
            &    1723805262919.0_dp/4571918432560.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,    &
            &    6338158500785.0_dp/11769362343261.0_dp, -4970555480458.0_dp/      &
            &    10924838743837.0_dp, 3326578051521.0_dp/2647936831840.0_dp,       &
            &    -880713585975.0_dp/1841400956686.0_dp, -1428733748635.0_dp/       &
            &    8843423958496.0_dp, 0.0_dp, 0.0_dp, 760814592956.0_dp/            &
            &    3276306540349.0_dp, 760814592956.0_dp/3276306540349.0_dp,         &
            &    -47223648122716.0_dp/6934462133451.0_dp, 71187472546993.0_dp/     &
            &    9669769126921.0_dp, -13330509492149.0_dp/9695768672337.0_dp,      &
            &    11565764226357.0_dp/8513123442827.0_dp, 0.0_dp], [7,7], order=[2,1])

            this%butcher_ass_exp(:)=this%butcher_ass_imp(:)
            this%butcher_c(:)=[0.247_dp, 4276536705230.0_dp/10142255878289.0_dp, &
            &                  0.335_dp, 0.075_dp, 0.7_dp, one]
