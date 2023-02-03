
@testset "$(rpad("Special Partitioned Additive Runge-Kutta Tableaus",80))" begin

    @test typeof(SPARKGLRK(1))                      <: SPARKMethod
    @test typeof(SPARKGLRK(2))                      <: SPARKMethod
    @test typeof(SPARKLobattoIIIAIIIB(2))           <: SPARKMethod
    @test typeof(SPARKLobattoIIIAIIIB(3))           <: SPARKMethod
    @test typeof(SPARKGLRKLobattoIIIAIIIB(1))       <: SPARKMethod
    @test typeof(SPARKGLRKLobattoIIIAIIIB(2))       <: SPARKMethod
    @test typeof(SPARKLobABC(2))                    <: SPARKMethod
    @test typeof(SPARKLobABC(3))                    <: SPARKMethod
    @test typeof(SPARKLobABD(2))                    <: SPARKMethod
    @test typeof(SPARKLobABD(3))                    <: SPARKMethod
    @test typeof(SPARKGLVPRK(1))                    <: SPARKMethod
    @test typeof(SPARKGLVPRK(2))                    <: SPARKMethod

    @test typeof(TableauGausspSymplectic(1))               <: VPARK
    @test typeof(TableauGausspSymplectic(2))               <: VPARK
    @test typeof(TableauLobattoIIIAIIIBpSymplectic(2))     <: VPARK
    @test typeof(TableauLobattoIIIAIIIBpSymplectic(3))     <: VPARK
    @test typeof(TableauLobattoIIIAIIIBpSymplectic(4))     <: VPARK
    @test typeof(TableauLobattoIIIBIIIApSymplectic(2))     <: VPARK
    @test typeof(TableauLobattoIIIBIIIApSymplectic(3))     <: VPARK
    @test typeof(TableauLobattoIIIBIIIApSymplectic(4))     <: VPARK

    @test typeof(TableauVSPARKGLRKpLobattoIIIAIIIB(1))                     <: VSPARKprimary
    @test typeof(TableauVSPARKGLRKpLobattoIIIAIIIB(2))                     <: VSPARKprimary
    @test typeof(TableauVSPARKGLRKpLobattoIIIBIIIA(1))                     <: VSPARKprimary
    @test typeof(TableauVSPARKGLRKpLobattoIIIBIIIA(2))                     <: VSPARKprimary
    @test typeof(TableauVSPARKGLRKpModifiedLobattoIIIAIIIB(1))             <: VSPARKprimary
    @test typeof(TableauVSPARKGLRKpModifiedLobattoIIIAIIIB(2))             <: VSPARKprimary
    @test typeof(TableauVSPARKGLRKpModifiedLobattoIIIBIIIA(1))             <: VSPARKprimary
    @test typeof(TableauVSPARKGLRKpModifiedLobattoIIIBIIIA(2))             <: VSPARKprimary
    @test typeof(TableauVSPARKGLRKpInternal(1))                            <: VSPARKprimary
    @test typeof(TableauVSPARKGLRKpInternal(2))                            <: VSPARKprimary
    @test typeof(TableauVSPARKGLRKpModifiedInternal(1))                    <: VSPARKprimary
    @test typeof(TableauVSPARKGLRKpModifiedInternal(2))                    <: VSPARKprimary
    @test typeof(TableauVSPARKGLRKpMidpoint(1))                            <: VSPARKprimary
    @test typeof(TableauVSPARKGLRKpMidpoint(2))                            <: VSPARKprimary
    @test typeof(TableauVSPARKGLRKpModifiedMidpoint(1))                    <: VSPARKprimary
    @test typeof(TableauVSPARKGLRKpModifiedMidpoint(2))                    <: VSPARKprimary
    @test typeof(TableauVSPARKGLRKpSymmetric(1))                           <: VSPARKprimary
    @test typeof(TableauVSPARKGLRKpSymmetric(2))                           <: VSPARKprimary
    @test typeof(TableauVSPARKGLRKpSymplectic(1))                          <: VSPARKprimary
    @test typeof(TableauVSPARKGLRKpSymplectic(2))                          <: VSPARKprimary

    @test typeof(TableauVSPARKLobattoIIIAIIIBpLobattoIIIAIIIB(2))          <: VSPARKprimary
    @test typeof(TableauVSPARKLobattoIIIAIIIBpLobattoIIIAIIIB(3))          <: VSPARKprimary
    @test typeof(TableauVSPARKLobattoIIIBIIIApLobattoIIIAIIIB(2))          <: VSPARKprimary
    @test typeof(TableauVSPARKLobattoIIIBIIIApLobattoIIIAIIIB(3))          <: VSPARKprimary
    @test typeof(TableauVSPARKLobattoIIIAIIIBpLobattoIIIBIIIA(2))          <: VSPARKprimary
    @test typeof(TableauVSPARKLobattoIIIAIIIBpLobattoIIIBIIIA(3))          <: VSPARKprimary
    @test typeof(TableauVSPARKLobattoIIIBIIIApLobattoIIIBIIIA(2))          <: VSPARKprimary
    @test typeof(TableauVSPARKLobattoIIIBIIIApLobattoIIIBIIIA(3))          <: VSPARKprimary

    @test typeof(TableauVSPARKLobattoIIIAIIIBpModifiedLobattoIIIAIIIB(2))  <: VSPARKprimary
    @test typeof(TableauVSPARKLobattoIIIAIIIBpModifiedLobattoIIIAIIIB(3))  <: VSPARKprimary
    @test typeof(TableauVSPARKLobattoIIIAIIIBpModifiedLobattoIIIBIIIA(2))  <: VSPARKprimary
    @test typeof(TableauVSPARKLobattoIIIAIIIBpModifiedLobattoIIIBIIIA(3))  <: VSPARKprimary
    @test typeof(TableauVSPARKLobattoIIIBIIIApModifiedLobattoIIIAIIIB(2))  <: VSPARKprimary
    @test typeof(TableauVSPARKLobattoIIIBIIIApModifiedLobattoIIIAIIIB(3))  <: VSPARKprimary
    @test typeof(TableauVSPARKLobattoIIIBIIIApModifiedLobattoIIIBIIIA(2))  <: VSPARKprimary
    @test typeof(TableauVSPARKLobattoIIIBIIIApModifiedLobattoIIIBIIIA(3))  <: VSPARKprimary

    @test typeof(TableauVSPARKLobattoIIIAIIIBpMidpoint(2))                 <: VSPARKprimary
    @test typeof(TableauVSPARKLobattoIIIAIIIBpMidpoint(3))                 <: VSPARKprimary
    @test typeof(TableauVSPARKLobattoIIIBIIIApMidpoint(2))                 <: VSPARKprimary
    @test typeof(TableauVSPARKLobattoIIIBIIIApMidpoint(3))                 <: VSPARKprimary
    @test typeof(TableauVSPARKLobattoIIIAIIIBpModifiedMidpoint(2))         <: VSPARKprimary
    @test typeof(TableauVSPARKLobattoIIIAIIIBpModifiedMidpoint(3))         <: VSPARKprimary
    @test typeof(TableauVSPARKLobattoIIIBIIIApModifiedMidpoint(2))         <: VSPARKprimary
    @test typeof(TableauVSPARKLobattoIIIBIIIApModifiedMidpoint(3))         <: VSPARKprimary
    @test typeof(TableauVSPARKLobattoIIIAIIIBpSymmetric(2))                <: VSPARKprimary
    @test typeof(TableauVSPARKLobattoIIIAIIIBpSymmetric(3))                <: VSPARKprimary
    @test typeof(TableauVSPARKLobattoIIIBIIIApSymmetric(2))                <: VSPARKprimary
    @test typeof(TableauVSPARKLobattoIIIBIIIApSymmetric(3))                <: VSPARKprimary

    @test typeof(TableauVSPARKLobABCCD(2))            <: VSPARKprimary
    @test typeof(TableauVSPARKLobABCCE(2))            <: VSPARKprimary
    @test typeof(TableauVSPARKLobABDE(2))             <: VSPARKprimary
    @test typeof(TableauVSPARKLobABED(2))             <: VSPARKprimary
    @test typeof(TableauVSPARKLobABD(2))              <: VSPARKprimary
    @test typeof(TableauVSPARKLobABE(2))              <: VSPARKprimary
    @test typeof(TableauVSPARKLobDE(2))               <: VSPARKprimary
    @test typeof(TableauVSPARKLobED(2))               <: VSPARKprimary

    @test typeof(TableauVSPARKLobattoIIIAB(2))        <: VSPARKsecondary
    @test typeof(TableauVSPARKLobattoIIIBA(2))        <: VSPARKsecondary
    @test typeof(TableauVSPARKLobattoIIICC̄(2))        <: VSPARKsecondary
    @test typeof(TableauVSPARKLobattoIIIC̄C(2))        <: VSPARKsecondary
    @test typeof(TableauVSPARKLobattoIIID(2))         <: VSPARKsecondary
    @test typeof(TableauVSPARKLobattoIIIE(2))         <: VSPARKsecondary
    @test typeof(TableauVSPARKGLRKLobattoIIIAB(2))    <: VSPARKsecondary
    @test typeof(TableauVSPARKGLRKLobattoIIIBA(2))    <: VSPARKsecondary
    @test typeof(TableauVSPARKGLRKLobattoIIICC̄(2))    <: VSPARKsecondary
    @test typeof(TableauVSPARKGLRKLobattoIIIC̄C(2))    <: VSPARKsecondary
    @test typeof(TableauVSPARKGLRKLobattoIIID(2))     <: VSPARKsecondary
    @test typeof(TableauVSPARKGLRKLobattoIIIE(2))     <: VSPARKsecondary

    @test typeof(TableauHPARKGLRK(1))                 <: HPARK
    @test typeof(TableauHPARKGLRK(2))                 <: HPARK
    @test typeof(TableauHPARKLobattoIIIAIIIB(2))      <: HPARK
    @test typeof(TableauHPARKLobattoIIIAIIIB(3))      <: HPARK
    @test typeof(TableauHPARKLobattoIIIBIIIA(2))      <: HPARK
    @test typeof(TableauHPARKLobattoIIIBIIIA(3))      <: HPARK

    # @test TableauHPARKGLRK(1) == getTableauHPARK(:HPARKGLRK1, TableauGauss(1), TableauGauss(1))
    # @test TableauHPARKGLRK(2) == getTableauHPARK(:HPARKGLRK2, TableauGauss(2), TableauGauss(2))

    @test typeof(TableauHSPARKGLRKpSymmetric(1))                       <: HSPARKprimary
    @test typeof(TableauHSPARKGLRKpSymmetric(2))                       <: HSPARKprimary
    @test typeof(TableauHSPARKLobattoIIIAIIIBpSymmetric(2))            <: HSPARKprimary
    @test typeof(TableauHSPARKLobattoIIIAIIIBpSymmetric(3))            <: HSPARKprimary
    @test typeof(TableauHSPARKLobattoIIIAIIIBpSymmetric(4))            <: HSPARKprimary
    @test typeof(TableauHSPARKLobattoIIIBIIIApSymmetric(2))            <: HSPARKprimary
    @test typeof(TableauHSPARKLobattoIIIBIIIApSymmetric(3))            <: HSPARKprimary
    @test typeof(TableauHSPARKLobattoIIIBIIIApSymmetric(4))            <: HSPARKprimary

    @test typeof(TableauHSPARKLobattoIIIAB(2))                         <: HSPARKsecondary
    @test typeof(TableauHSPARKLobattoIIIBA(2))                         <: HSPARKsecondary
    @test typeof(TableauHSPARKLobattoIIICC̄(2))                         <: HSPARKsecondary
    @test typeof(TableauHSPARKLobattoIIIC̄C(2))                         <: HSPARKsecondary
    @test typeof(TableauHSPARKLobattoIIID(2))                          <: HSPARKsecondary
    @test typeof(TableauHSPARKLobattoIIIE(2))                          <: HSPARKsecondary
    @test typeof(TableauHSPARKGLRKLobattoIIIAB(1))                     <: HSPARKsecondary
    @test typeof(TableauHSPARKGLRKLobattoIIIBA(1))                     <: HSPARKsecondary
    @test typeof(TableauHSPARKGLRKLobattoIIICC̄(1))                     <: HSPARKsecondary
    @test typeof(TableauHSPARKGLRKLobattoIIIC̄C(1))                     <: HSPARKsecondary
    @test typeof(TableauHSPARKGLRKLobattoIIID(1))                      <: HSPARKsecondary
    @test typeof(TableauHSPARKGLRKLobattoIIIE(1))                      <: HSPARKsecondary

    @test typeof(SLRKLobattoIIIAB(2))          <: SLRK
    @test typeof(SLRKLobattoIIIBA(2))          <: SLRK
    @test typeof(SLRKLobattoIIICC̄(2))          <: SLRK
    @test typeof(SLRKLobattoIIIC̄C(2))          <: SLRK
    @test typeof(SLRKLobattoIIID(2))           <: SLRK
    @test typeof(SLRKLobattoIIIE(2))           <: SLRK

end
