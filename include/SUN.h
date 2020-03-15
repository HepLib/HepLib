Tensor T,f(antisymmetric);
Tensor colTp;
Symbols I2R,NF,NA;
Dimension NA;
AutoDeclare Index colAj;
Dimension NF;
AutoDeclare Index colFi;

#procedure SUNTrace

#do colXi = 1,1
    if ( count(f,1) || match(T(colFi1?,colFi2?,colAj1?)*T(colFi3?,colFi4?,colAj1?)) ) redefine colXi "0";
    id,once,f(colAj1?,colAj2?,colAj3?) = 1/I2R/i_*T(colFi1,colFi2,colAj1)*T(colFi2,colFi3,colAj2)*T(colFi3,colFi1,colAj3)-1/I2R/i_*T(colFi1,colFi2,colAj3)*T(colFi2,colFi3,colAj2)*T(colFi3,colFi1,colAj1);
    sum colFi1,colFi2,colFi3;
    id T(colFi1?,colFi2?,colAj1?)*T(colFi3?,colFi4?,colAj1?) = colTp(colFi1,colFi2,colFi3,colFi4);
    #do colXj = 1,1
        if ( count(colTp,1) ) redefine colXj "0";
        .sort
        id,once,colTp(colFi1?,colFi2?,colFi3?,colFi4?) = I2R*(d_(colFi1,colFi4)*d_(colFi2,colFi3)-d_(colFi1,colFi2)*d_(colFi3,colFi4)/NF);
    #enddo
#enddo

#endprocedure
