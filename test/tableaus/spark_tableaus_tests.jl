
@testset "$(rpad("Special Partitioned Additive Runge-Kutta Tableaus",80))" begin

    @test typeof(TableauSPARKGLRK(1))                      <: TableauSPARK
    @test typeof(TableauSPARKGLRK(2))                      <: TableauSPARK
    @test typeof(TableauSPARKLobattoIIIAIIIB(2))           <: TableauSPARK
    @test typeof(TableauSPARKLobattoIIIAIIIB(3))           <: TableauSPARK
    @test typeof(TableauSPARKGLRKLobattoIIIAIIIB(1))       <: TableauSPARK
    @test typeof(TableauSPARKGLRKLobattoIIIAIIIB(2))       <: TableauSPARK
    @test typeof(TableauSPARKLobABC(2))                    <: TableauSPARK
    @test typeof(TableauSPARKLobABC(3))                    <: TableauSPARK
    @test typeof(TableauSPARKLobABD(2))                    <: TableauSPARK
    @test typeof(TableauSPARKLobABD(3))                    <: TableauSPARK
    @test typeof(TableauSPARKGLVPRK(1))                    <: TableauSPARK
    @test typeof(TableauSPARKGLVPRK(2))                    <: TableauSPARK

    @test typeof(TableauGLRKpSymplectic(1))                <: TableauVPARK
    @test typeof(TableauGLRKpSymplectic(2))                <: TableauVPARK
    @test typeof(TableauLobattoIIIAIIIBpSymplectic(2))     <: TableauVPARK
    @test typeof(TableauLobattoIIIAIIIBpSymplectic(3))     <: TableauVPARK
    @test typeof(TableauLobattoIIIAIIIBpSymplectic(4))     <: TableauVPARK
    @test typeof(TableauLobattoIIIBIIIApSymplectic(2))     <: TableauVPARK
    @test typeof(TableauLobattoIIIBIIIApSymplectic(3))     <: TableauVPARK
    @test typeof(TableauLobattoIIIBIIIApSymplectic(4))     <: TableauVPARK

    @test typeof(TableauVSPARKGLRKpLobattoIIIAIIIB(1))                     <: TableauVSPARKprimary
    @test typeof(TableauVSPARKGLRKpLobattoIIIAIIIB(2))                     <: TableauVSPARKprimary
    @test typeof(TableauVSPARKGLRKpLobattoIIIBIIIA(1))                     <: TableauVSPARKprimary
    @test typeof(TableauVSPARKGLRKpLobattoIIIBIIIA(2))                     <: TableauVSPARKprimary
    @test typeof(TableauVSPARKGLRKpModifiedLobattoIIIAIIIB(1))             <: TableauVSPARKprimary
    @test typeof(TableauVSPARKGLRKpModifiedLobattoIIIAIIIB(2))             <: TableauVSPARKprimary
    @test typeof(TableauVSPARKGLRKpModifiedLobattoIIIBIIIA(1))             <: TableauVSPARKprimary
    @test typeof(TableauVSPARKGLRKpModifiedLobattoIIIBIIIA(2))             <: TableauVSPARKprimary
    @test typeof(TableauVSPARKGLRKpInternal(1))                            <: TableauVSPARKprimary
    @test typeof(TableauVSPARKGLRKpInternal(2))                            <: TableauVSPARKprimary
    @test typeof(TableauVSPARKGLRKpModifiedInternal(1))                    <: TableauVSPARKprimary
    @test typeof(TableauVSPARKGLRKpModifiedInternal(2))                    <: TableauVSPARKprimary
    @test typeof(TableauVSPARKGLRKpMidpoint(1))                            <: TableauVSPARKprimary
    @test typeof(TableauVSPARKGLRKpMidpoint(2))                            <: TableauVSPARKprimary
    @test typeof(TableauVSPARKGLRKpModifiedMidpoint(1))                    <: TableauVSPARKprimary
    @test typeof(TableauVSPARKGLRKpModifiedMidpoint(2))                    <: TableauVSPARKprimary
    @test typeof(TableauVSPARKGLRKpSymmetric(1))                           <: TableauVSPARKprimary
    @test typeof(TableauVSPARKGLRKpSymmetric(2))                           <: TableauVSPARKprimary
    @test typeof(TableauVSPARKGLRKpSymplectic(1))                          <: TableauVSPARKprimary
    @test typeof(TableauVSPARKGLRKpSymplectic(2))                          <: TableauVSPARKprimary

    @test typeof(TableauVSPARKLobattoIIIAIIIBpLobattoIIIAIIIB(2))          <: TableauVSPARKprimary
    @test typeof(TableauVSPARKLobattoIIIAIIIBpLobattoIIIAIIIB(3))          <: TableauVSPARKprimary
    @test typeof(TableauVSPARKLobattoIIIBIIIApLobattoIIIAIIIB(2))          <: TableauVSPARKprimary
    @test typeof(TableauVSPARKLobattoIIIBIIIApLobattoIIIAIIIB(3))          <: TableauVSPARKprimary
    @test typeof(TableauVSPARKLobattoIIIAIIIBpLobattoIIIBIIIA(2))          <: TableauVSPARKprimary
    @test typeof(TableauVSPARKLobattoIIIAIIIBpLobattoIIIBIIIA(3))          <: TableauVSPARKprimary
    @test typeof(TableauVSPARKLobattoIIIBIIIApLobattoIIIBIIIA(2))          <: TableauVSPARKprimary
    @test typeof(TableauVSPARKLobattoIIIBIIIApLobattoIIIBIIIA(3))          <: TableauVSPARKprimary

    @test typeof(TableauVSPARKLobattoIIIAIIIBpModifiedLobattoIIIAIIIB(2))  <: TableauVSPARKprimary
    @test typeof(TableauVSPARKLobattoIIIAIIIBpModifiedLobattoIIIAIIIB(3))  <: TableauVSPARKprimary
    @test typeof(TableauVSPARKLobattoIIIAIIIBpModifiedLobattoIIIBIIIA(2))  <: TableauVSPARKprimary
    @test typeof(TableauVSPARKLobattoIIIAIIIBpModifiedLobattoIIIBIIIA(3))  <: TableauVSPARKprimary
    @test typeof(TableauVSPARKLobattoIIIBIIIApModifiedLobattoIIIAIIIB(2))  <: TableauVSPARKprimary
    @test typeof(TableauVSPARKLobattoIIIBIIIApModifiedLobattoIIIAIIIB(3))  <: TableauVSPARKprimary
    @test typeof(TableauVSPARKLobattoIIIBIIIApModifiedLobattoIIIBIIIA(2))  <: TableauVSPARKprimary
    @test typeof(TableauVSPARKLobattoIIIBIIIApModifiedLobattoIIIBIIIA(3))  <: TableauVSPARKprimary

    @test typeof(TableauVSPARKLobattoIIIAIIIBpMidpoint(2))                 <: TableauVSPARKprimary
    @test typeof(TableauVSPARKLobattoIIIAIIIBpMidpoint(3))                 <: TableauVSPARKprimary
    @test typeof(TableauVSPARKLobattoIIIBIIIApMidpoint(2))                 <: TableauVSPARKprimary
    @test typeof(TableauVSPARKLobattoIIIBIIIApMidpoint(3))                 <: TableauVSPARKprimary
    @test typeof(TableauVSPARKLobattoIIIAIIIBpModifiedMidpoint(2))         <: TableauVSPARKprimary
    @test typeof(TableauVSPARKLobattoIIIAIIIBpModifiedMidpoint(3))         <: TableauVSPARKprimary
    @test typeof(TableauVSPARKLobattoIIIBIIIApModifiedMidpoint(2))         <: TableauVSPARKprimary
    @test typeof(TableauVSPARKLobattoIIIBIIIApModifiedMidpoint(3))         <: TableauVSPARKprimary
    @test typeof(TableauVSPARKLobattoIIIAIIIBpSymmetric(2))                <: TableauVSPARKprimary
    @test typeof(TableauVSPARKLobattoIIIAIIIBpSymmetric(3))                <: TableauVSPARKprimary
    @test typeof(TableauVSPARKLobattoIIIBIIIApSymmetric(2))                <: TableauVSPARKprimary
    @test typeof(TableauVSPARKLobattoIIIBIIIApSymmetric(3))                <: TableauVSPARKprimary

    @test typeof(TableauVSPARKLobABCCD(2))            <: TableauVSPARKprimary
    @test typeof(TableauVSPARKLobABCCE(2))            <: TableauVSPARKprimary
    @test typeof(TableauVSPARKLobABDE(2))             <: TableauVSPARKprimary
    @test typeof(TableauVSPARKLobABED(2))             <: TableauVSPARKprimary
    @test typeof(TableauVSPARKLobABD(2))              <: TableauVSPARKprimary
    @test typeof(TableauVSPARKLobABE(2))              <: TableauVSPARKprimary
    @test typeof(TableauVSPARKLobDE(2))               <: TableauVSPARKprimary
    @test typeof(TableauVSPARKLobED(2))               <: TableauVSPARKprimary

    @test typeof(TableauVSPARKLobattoIIIAB(2))        <: TableauVSPARKsecondary
    @test typeof(TableauVSPARKLobattoIIIBA(2))        <: TableauVSPARKsecondary
    @test typeof(TableauVSPARKLobattoIIICC̄(2))        <: TableauVSPARKsecondary
    @test typeof(TableauVSPARKLobattoIIIC̄C(2))        <: TableauVSPARKsecondary
    @test typeof(TableauVSPARKLobattoIIID(2))         <: TableauVSPARKsecondary
    @test typeof(TableauVSPARKLobattoIIIE(2))         <: TableauVSPARKsecondary
    @test typeof(TableauVSPARKGLRKLobattoIIIAB(2))    <: TableauVSPARKsecondary
    @test typeof(TableauVSPARKGLRKLobattoIIIBA(2))    <: TableauVSPARKsecondary
    @test typeof(TableauVSPARKGLRKLobattoIIICC̄(2))    <: TableauVSPARKsecondary
    @test typeof(TableauVSPARKGLRKLobattoIIIC̄C(2))    <: TableauVSPARKsecondary
    @test typeof(TableauVSPARKGLRKLobattoIIID(2))     <: TableauVSPARKsecondary
    @test typeof(TableauVSPARKGLRKLobattoIIIE(2))     <: TableauVSPARKsecondary

    @test typeof(TableauHPARKGLRK(1))                 <: TableauHPARK
    @test typeof(TableauHPARKGLRK(2))                 <: TableauHPARK
    @test typeof(TableauHPARKLobattoIIIAIIIB(2))      <: TableauHPARK
    @test typeof(TableauHPARKLobattoIIIAIIIB(3))      <: TableauHPARK
    @test typeof(TableauHPARKLobattoIIIBIIIA(2))      <: TableauHPARK
    @test typeof(TableauHPARKLobattoIIIBIIIA(3))      <: TableauHPARK

    @test TableauHPARKGLRK(1) == getTableauHPARK(:HPARKGLRK1, CoefficientsGLRK(1), CoefficientsGLRK(1))
    @test TableauHPARKGLRK(2) == getTableauHPARK(:HPARKGLRK2, CoefficientsGLRK(2), CoefficientsGLRK(2))

    @test typeof(TableauHSPARKGLRKpSymmetric(1))                       <: TableauHSPARKprimary
    @test typeof(TableauHSPARKGLRKpSymmetric(2))                       <: TableauHSPARKprimary
    @test typeof(TableauHSPARKLobattoIIIAIIIBpSymmetric(2))            <: TableauHSPARKprimary
    @test typeof(TableauHSPARKLobattoIIIAIIIBpSymmetric(3))            <: TableauHSPARKprimary
    @test typeof(TableauHSPARKLobattoIIIAIIIBpSymmetric(4))            <: TableauHSPARKprimary
    @test typeof(TableauHSPARKLobattoIIIBIIIApSymmetric(2))            <: TableauHSPARKprimary
    @test typeof(TableauHSPARKLobattoIIIBIIIApSymmetric(3))            <: TableauHSPARKprimary
    @test typeof(TableauHSPARKLobattoIIIBIIIApSymmetric(4))            <: TableauHSPARKprimary

    @test typeof(TableauHSPARKLobattoIIIAB(2))                         <: TableauHSPARKsecondary
    @test typeof(TableauHSPARKLobattoIIIBA(2))                         <: TableauHSPARKsecondary
    @test typeof(TableauHSPARKLobattoIIICC̄(2))                         <: TableauHSPARKsecondary
    @test typeof(TableauHSPARKLobattoIIIC̄C(2))                         <: TableauHSPARKsecondary
    @test typeof(TableauHSPARKLobattoIIID(2))                          <: TableauHSPARKsecondary
    @test typeof(TableauHSPARKLobattoIIIE(2))                          <: TableauHSPARKsecondary
    @test typeof(TableauHSPARKGLRKLobattoIIIAB(1))                     <: TableauHSPARKsecondary
    @test typeof(TableauHSPARKGLRKLobattoIIIBA(1))                     <: TableauHSPARKsecondary
    @test typeof(TableauHSPARKGLRKLobattoIIICC̄(1))                     <: TableauHSPARKsecondary
    @test typeof(TableauHSPARKGLRKLobattoIIIC̄C(1))                     <: TableauHSPARKsecondary
    @test typeof(TableauHSPARKGLRKLobattoIIID(1))                      <: TableauHSPARKsecondary
    @test typeof(TableauHSPARKGLRKLobattoIIIE(1))                      <: TableauHSPARKsecondary

    @test typeof(TableauSLRKLobattoIIIAB(2))          <: TableauSLRK
    @test typeof(TableauSLRKLobattoIIIBA(2))          <: TableauSLRK
    @test typeof(TableauSLRKLobattoIIICC̄(2))          <: TableauSLRK
    @test typeof(TableauSLRKLobattoIIIC̄C(2))          <: TableauSLRK
    @test typeof(TableauSLRKLobattoIIID(2))           <: TableauSLRK
    @test typeof(TableauSLRKLobattoIIIE(2))           <: TableauSLRK

end
