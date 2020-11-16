
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

    @test typeof(getTableauVSPARKGLRKpLobIIIAIIIB(1))             <: TableauVSPARKprimary
    @test typeof(getTableauVSPARKGLRKpLobIIIAIIIB(2))             <: TableauVSPARKprimary
    @test typeof(getTableauVSPARKGLRKpLobIIIBIIIA(1))             <: TableauVSPARKprimary
    @test typeof(getTableauVSPARKGLRKpLobIIIBIIIA(2))             <: TableauVSPARKprimary
    @test typeof(getTableauVSPARKGLRKpModifiedLobIIIAIIIB(1))     <: TableauVSPARKprimary
    @test typeof(getTableauVSPARKGLRKpModifiedLobIIIAIIIB(2))     <: TableauVSPARKprimary
    @test typeof(getTableauVSPARKGLRKpModifiedLobIIIBIIIA(1))     <: TableauVSPARKprimary
    @test typeof(getTableauVSPARKGLRKpModifiedLobIIIBIIIA(2))     <: TableauVSPARKprimary
    @test typeof(getTableauVSPARKGLRKpInternal(1))                <: TableauVSPARKprimary
    @test typeof(getTableauVSPARKGLRKpInternal(2))                <: TableauVSPARKprimary
    @test typeof(getTableauVSPARKGLRKpModifiedInternal(1))        <: TableauVSPARKprimary
    @test typeof(getTableauVSPARKGLRKpModifiedInternal(2))        <: TableauVSPARKprimary
    @test typeof(getTableauVSPARKGLRKpMidpoint(1))                <: TableauVSPARKprimary
    @test typeof(getTableauVSPARKGLRKpMidpoint(2))                <: TableauVSPARKprimary
    @test typeof(getTableauVSPARKGLRKpModifiedMidpoint(1))        <: TableauVSPARKprimary
    @test typeof(getTableauVSPARKGLRKpModifiedMidpoint(2))        <: TableauVSPARKprimary
    @test typeof(getTableauVSPARKGLRKpSymmetric(1))               <: TableauVSPARKprimary
    @test typeof(getTableauVSPARKGLRKpSymmetric(2))               <: TableauVSPARKprimary
    @test typeof(getTableauVSPARKGLRKpSymplectic(1))              <: TableauVSPARKprimary
    @test typeof(getTableauVSPARKGLRKpSymplectic(2))              <: TableauVSPARKprimary

    @test typeof(getTableauVSPARKLobIIIAIIIBpLobIIIAIIIB(2))      <: TableauVSPARKprimary
    @test typeof(getTableauVSPARKLobIIIAIIIBpLobIIIAIIIB(3))      <: TableauVSPARKprimary
    @test typeof(getTableauVSPARKLobIIIBIIIApLobIIIAIIIB(2))      <: TableauVSPARKprimary
    @test typeof(getTableauVSPARKLobIIIBIIIApLobIIIAIIIB(3))      <: TableauVSPARKprimary
    @test typeof(getTableauVSPARKLobIIIAIIIBpLobIIIBIIIA(2))      <: TableauVSPARKprimary
    @test typeof(getTableauVSPARKLobIIIAIIIBpLobIIIBIIIA(3))      <: TableauVSPARKprimary
    @test typeof(getTableauVSPARKLobIIIBIIIApLobIIIBIIIA(2))      <: TableauVSPARKprimary
    @test typeof(getTableauVSPARKLobIIIBIIIApLobIIIBIIIA(3))      <: TableauVSPARKprimary

    @test typeof(getTableauVSPARKLobIIIAIIIBpModifiedLobIIIAIIIB(2))      <: TableauVSPARKprimary
    @test typeof(getTableauVSPARKLobIIIAIIIBpModifiedLobIIIAIIIB(3))      <: TableauVSPARKprimary
    @test typeof(getTableauVSPARKLobIIIAIIIBpModifiedLobIIIBIIIA(2))      <: TableauVSPARKprimary
    @test typeof(getTableauVSPARKLobIIIAIIIBpModifiedLobIIIBIIIA(3))      <: TableauVSPARKprimary
    @test typeof(getTableauVSPARKLobIIIBIIIApModifiedLobIIIAIIIB(2))      <: TableauVSPARKprimary
    @test typeof(getTableauVSPARKLobIIIBIIIApModifiedLobIIIAIIIB(3))      <: TableauVSPARKprimary
    @test typeof(getTableauVSPARKLobIIIBIIIApModifiedLobIIIBIIIA(2))      <: TableauVSPARKprimary
    @test typeof(getTableauVSPARKLobIIIBIIIApModifiedLobIIIBIIIA(3))      <: TableauVSPARKprimary

    @test typeof(getTableauVSPARKLobIIIAIIIBpMidpoint(2))         <: TableauVSPARKprimary
    @test typeof(getTableauVSPARKLobIIIAIIIBpMidpoint(3))         <: TableauVSPARKprimary
    @test typeof(getTableauVSPARKLobIIIBIIIApMidpoint(2))         <: TableauVSPARKprimary
    @test typeof(getTableauVSPARKLobIIIBIIIApMidpoint(3))         <: TableauVSPARKprimary
    @test typeof(getTableauVSPARKLobIIIAIIIBpModifiedMidpoint(2)) <: TableauVSPARKprimary
    @test typeof(getTableauVSPARKLobIIIAIIIBpModifiedMidpoint(3)) <: TableauVSPARKprimary
    @test typeof(getTableauVSPARKLobIIIBIIIApModifiedMidpoint(2)) <: TableauVSPARKprimary
    @test typeof(getTableauVSPARKLobIIIBIIIApModifiedMidpoint(3)) <: TableauVSPARKprimary
    @test typeof(getTableauVSPARKLobIIIAIIIBpSymmetric(2))        <: TableauVSPARKprimary
    @test typeof(getTableauVSPARKLobIIIAIIIBpSymmetric(3))        <: TableauVSPARKprimary
    @test typeof(getTableauVSPARKLobIIIBIIIApSymmetric(2))        <: TableauVSPARKprimary
    @test typeof(getTableauVSPARKLobIIIBIIIApSymmetric(3))        <: TableauVSPARKprimary

    @test typeof(getTableauVSPARKLobABCCD(2))     <: TableauVSPARKprimary
    @test typeof(getTableauVSPARKLobABCCE(2))     <: TableauVSPARKprimary
    @test typeof(getTableauVSPARKLobABDE(2))      <: TableauVSPARKprimary
    @test typeof(getTableauVSPARKLobABED(2))      <: TableauVSPARKprimary
    @test typeof(getTableauVSPARKLobABD(2))       <: TableauVSPARKprimary
    @test typeof(getTableauVSPARKLobABE(2))       <: TableauVSPARKprimary
    @test typeof(getTableauVSPARKLobDE(2))        <: TableauVSPARKprimary
    @test typeof(getTableauVSPARKLobED(2))        <: TableauVSPARKprimary

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
