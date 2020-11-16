
@testset "$(rpad("Special Partitioned Additive Runge-Kutta Tableaus",80))" begin

    @test typeof(TableauSPARKGLRK(1))                 <: TableauSPARK
    @test typeof(TableauSPARKGLRK(2))                 <: TableauSPARK
    @test typeof(TableauSPARKLobIIIAIIIB(2))          <: TableauSPARK
    @test typeof(TableauSPARKLobIIIAIIIB(3))          <: TableauSPARK
    @test typeof(TableauSPARKGLRKLobIIIAIIIB(1))      <: TableauSPARK
    @test typeof(TableauSPARKGLRKLobIIIAIIIB(2))      <: TableauSPARK
    @test typeof(TableauSPARKLobABC(2))               <: TableauSPARK
    @test typeof(TableauSPARKLobABC(3))               <: TableauSPARK
    @test typeof(TableauSPARKLobABD(2))               <: TableauSPARK
    @test typeof(TableauSPARKLobABD(3))               <: TableauSPARK
    @test typeof(TableauSPARKGLVPRK(1))               <: TableauSPARK
    @test typeof(TableauSPARKGLVPRK(2))               <: TableauSPARK

    @test typeof(getTableauGLRKpSymplectic(1))        <: TableauVPARK
    @test typeof(getTableauGLRKpSymplectic(2))        <: TableauVPARK
    @test typeof(getTableauLobIIIAIIIB2pSymplectic()) <: TableauVPARK
    @test typeof(getTableauLobIIIAIIIB3pSymplectic()) <: TableauVPARK
    @test typeof(getTableauLobIIIAIIIB4pSymplectic()) <: TableauVPARK

    @test typeof(TableauVSPARKGLRKpLobIIIAIIIB(1))                 <: TableauVSPARKprimary
    @test typeof(TableauVSPARKGLRKpLobIIIAIIIB(2))                 <: TableauVSPARKprimary
    @test typeof(TableauVSPARKGLRKpLobIIIBIIIA(1))                 <: TableauVSPARKprimary
    @test typeof(TableauVSPARKGLRKpLobIIIBIIIA(2))                 <: TableauVSPARKprimary
    @test typeof(TableauVSPARKGLRKpModifiedLobIIIAIIIB(1))         <: TableauVSPARKprimary
    @test typeof(TableauVSPARKGLRKpModifiedLobIIIAIIIB(2))         <: TableauVSPARKprimary
    @test typeof(TableauVSPARKGLRKpModifiedLobIIIBIIIA(1))         <: TableauVSPARKprimary
    @test typeof(TableauVSPARKGLRKpModifiedLobIIIBIIIA(2))         <: TableauVSPARKprimary
    @test typeof(TableauVSPARKGLRKpInternal(1))                    <: TableauVSPARKprimary
    @test typeof(TableauVSPARKGLRKpInternal(2))                    <: TableauVSPARKprimary
    @test typeof(TableauVSPARKGLRKpModifiedInternal(1))            <: TableauVSPARKprimary
    @test typeof(TableauVSPARKGLRKpModifiedInternal(2))            <: TableauVSPARKprimary
    @test typeof(TableauVSPARKGLRKpMidpoint(1))                    <: TableauVSPARKprimary
    @test typeof(TableauVSPARKGLRKpMidpoint(2))                    <: TableauVSPARKprimary
    @test typeof(TableauVSPARKGLRKpModifiedMidpoint(1))            <: TableauVSPARKprimary
    @test typeof(TableauVSPARKGLRKpModifiedMidpoint(2))            <: TableauVSPARKprimary
    @test typeof(TableauVSPARKGLRKpSymmetric(1))                   <: TableauVSPARKprimary
    @test typeof(TableauVSPARKGLRKpSymmetric(2))                   <: TableauVSPARKprimary
    @test typeof(TableauVSPARKGLRKpSymplectic(1))                  <: TableauVSPARKprimary
    @test typeof(TableauVSPARKGLRKpSymplectic(2))                  <: TableauVSPARKprimary

    @test typeof(TableauVSPARKLobIIIAIIIBpLobIIIAIIIB(2))          <: TableauVSPARKprimary
    @test typeof(TableauVSPARKLobIIIAIIIBpLobIIIAIIIB(3))          <: TableauVSPARKprimary
    @test typeof(TableauVSPARKLobIIIBIIIApLobIIIAIIIB(2))          <: TableauVSPARKprimary
    @test typeof(TableauVSPARKLobIIIBIIIApLobIIIAIIIB(3))          <: TableauVSPARKprimary
    @test typeof(TableauVSPARKLobIIIAIIIBpLobIIIBIIIA(2))          <: TableauVSPARKprimary
    @test typeof(TableauVSPARKLobIIIAIIIBpLobIIIBIIIA(3))          <: TableauVSPARKprimary
    @test typeof(TableauVSPARKLobIIIBIIIApLobIIIBIIIA(2))          <: TableauVSPARKprimary
    @test typeof(TableauVSPARKLobIIIBIIIApLobIIIBIIIA(3))          <: TableauVSPARKprimary

    @test typeof(TableauVSPARKLobIIIAIIIBpModifiedLobIIIAIIIB(2))  <: TableauVSPARKprimary
    @test typeof(TableauVSPARKLobIIIAIIIBpModifiedLobIIIAIIIB(3))  <: TableauVSPARKprimary
    @test typeof(TableauVSPARKLobIIIAIIIBpModifiedLobIIIBIIIA(2))  <: TableauVSPARKprimary
    @test typeof(TableauVSPARKLobIIIAIIIBpModifiedLobIIIBIIIA(3))  <: TableauVSPARKprimary
    @test typeof(TableauVSPARKLobIIIBIIIApModifiedLobIIIAIIIB(2))  <: TableauVSPARKprimary
    @test typeof(TableauVSPARKLobIIIBIIIApModifiedLobIIIAIIIB(3))  <: TableauVSPARKprimary
    @test typeof(TableauVSPARKLobIIIBIIIApModifiedLobIIIBIIIA(2))  <: TableauVSPARKprimary
    @test typeof(TableauVSPARKLobIIIBIIIApModifiedLobIIIBIIIA(3))  <: TableauVSPARKprimary

    @test typeof(TableauVSPARKLobIIIAIIIBpMidpoint(2))             <: TableauVSPARKprimary
    @test typeof(TableauVSPARKLobIIIAIIIBpMidpoint(3))             <: TableauVSPARKprimary
    @test typeof(TableauVSPARKLobIIIBIIIApMidpoint(2))             <: TableauVSPARKprimary
    @test typeof(TableauVSPARKLobIIIBIIIApMidpoint(3))             <: TableauVSPARKprimary
    @test typeof(TableauVSPARKLobIIIAIIIBpModifiedMidpoint(2))     <: TableauVSPARKprimary
    @test typeof(TableauVSPARKLobIIIAIIIBpModifiedMidpoint(3))     <: TableauVSPARKprimary
    @test typeof(TableauVSPARKLobIIIBIIIApModifiedMidpoint(2))     <: TableauVSPARKprimary
    @test typeof(TableauVSPARKLobIIIBIIIApModifiedMidpoint(3))     <: TableauVSPARKprimary
    @test typeof(TableauVSPARKLobIIIAIIIBpSymmetric(2))            <: TableauVSPARKprimary
    @test typeof(TableauVSPARKLobIIIAIIIBpSymmetric(3))            <: TableauVSPARKprimary
    @test typeof(TableauVSPARKLobIIIBIIIApSymmetric(2))            <: TableauVSPARKprimary
    @test typeof(TableauVSPARKLobIIIBIIIApSymmetric(3))            <: TableauVSPARKprimary

    @test typeof(TableauVSPARKLobABCCD(2))        <: TableauVSPARKprimary
    @test typeof(TableauVSPARKLobABCCE(2))        <: TableauVSPARKprimary
    @test typeof(TableauVSPARKLobABDE(2))         <: TableauVSPARKprimary
    @test typeof(TableauVSPARKLobABED(2))         <: TableauVSPARKprimary
    @test typeof(TableauVSPARKLobABD(2))          <: TableauVSPARKprimary
    @test typeof(TableauVSPARKLobABE(2))          <: TableauVSPARKprimary
    @test typeof(TableauVSPARKLobDE(2))           <: TableauVSPARKprimary
    @test typeof(TableauVSPARKLobED(2))           <: TableauVSPARKprimary

    @test typeof(getTableauVSPARKLobIIIAB(2))     <: TableauVSPARKsecondary
    @test typeof(getTableauVSPARKLobIIIC(2))      <: TableauVSPARKsecondary
    @test typeof(getTableauVSPARKLobIIID(2))      <: TableauVSPARKsecondary
    @test typeof(getTableauVSPARKLobIIIE(2))      <: TableauVSPARKsecondary
    @test typeof(getTableauVSPARKGLRKLobIIIAB(2)) <: TableauVSPARKsecondary
    @test typeof(getTableauVSPARKGLRKLobIIIC(2))  <: TableauVSPARKsecondary
    @test typeof(getTableauVSPARKGLRKLobIIID(2))  <: TableauVSPARKsecondary
    @test typeof(getTableauVSPARKGLRKLobIIIE(2))  <: TableauVSPARKsecondary

    @test typeof(getTableauHPARKGLRK(1))          <: TableauHPARK
    @test typeof(getTableauHPARKGLRK(2))          <: TableauHPARK
    @test typeof(getTableauHPARKLobIIIAIIIB2())   <: TableauHPARK
    @test typeof(getTableauHPARKLobIIIAIIIB3())   <: TableauHPARK
    @test typeof(getTableauHPARKLobIIIAIIIB4())   <: TableauHPARK

    @test typeof(getTableauHSPARKGLRKpSymmetric(1))        <: TableauHSPARKprimary
    @test typeof(getTableauHSPARKGLRKpSymmetric(2))        <: TableauHSPARKprimary
    @test typeof(getTableauHSPARKLobIIIAIIIB2pSymmetric()) <: TableauHSPARKprimary
    @test typeof(getTableauHSPARKLobIIIAIIIB3pSymmetric()) <: TableauHSPARKprimary
    @test typeof(getTableauHSPARKLobIIIAIIIB4pSymmetric()) <: TableauHSPARKprimary

    @test typeof(getTableauHSPARKLobIIIAB(2))     <: TableauHSPARKsecondary
    @test typeof(getTableauHSPARKLobIIIC(2))      <: TableauHSPARKsecondary
    @test typeof(getTableauHSPARKLobIIID(2))      <: TableauHSPARKsecondary
    @test typeof(getTableauHSPARKLobIIIE(2))      <: TableauHSPARKsecondary
    @test typeof(getTableauHSPARKGLRKLobIIIAB(1)) <: TableauHSPARKsecondary
    @test typeof(getTableauHSPARKGLRKLobIIIC(1))  <: TableauHSPARKsecondary
    @test typeof(getTableauHSPARKGLRKLobIIID(1))  <: TableauHSPARKsecondary
    @test typeof(getTableauHSPARKGLRKLobIIIE(1))  <: TableauHSPARKsecondary

    @test typeof(getTableauSLRKLobIIIAB(2)) <: TableauSLRK
    @test typeof(getTableauSLRKLobIIIC(2))  <: TableauSLRK
    @test typeof(getTableauSLRKLobIIID(2))  <: TableauSLRK
    @test typeof(getTableauSLRKLobIIIE(2))  <: TableauSLRK

end
