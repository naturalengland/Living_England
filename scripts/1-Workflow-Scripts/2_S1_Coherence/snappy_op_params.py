"""-------------------------------------------------------------
    Script Name:   snappy_op_params.py
    Description:   Identifies SNAP operator parameters and default values
    Created By:    Chris Moore
    Date Created:  10.12.21
    Date Modified: 26.07.22
-------------------------------------------------------------"""

import sys

print('> snappy_op_params\n')
# input snap operator
op = sys.argv[1]
# catch format error
if op[0] == "'":
    print("ERROR: '' not supported")

# load snap and run
else:
    from snappy import GPF

    # get operator
    op_spi = GPF.getDefaultInstance().getOperatorSpiRegistry().getOperatorSpi(op)
    # catch operator name error
    if not op_spi:
        print('\nERROR: invalid operator')
    else:
        # get descriptors
        param_Desc = op_spi.getOperatorDescriptor().getParameterDescriptors()
        # print operator and descriptor info
        print('\n> Operator name: {}\n> Operator alias: {}\n'.format(op_spi.getOperatorDescriptor().getName(),
                                                                     op_spi.getOperatorDescriptor().getAlias()))
        for param in param_Desc:
            print('> Parameter name: {}\n> Parameter alias: {}\n> Default value: {}\n'.format(param.getName(),
                                                                                              param.getAlias(),
                                                                                              param.getDefaultValue()))
