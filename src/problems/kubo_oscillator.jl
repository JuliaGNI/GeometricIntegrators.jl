module KuboOscillatorProblem

    using GeometricIntegrators.Equations

    export kubo_oscillator_sde_1, kubo_oscillator_psde_1, kubo_oscillator_spsde_1
    export kubo_oscillator_sde_2, kubo_oscillator_psde_2, kubo_oscillator_spsde_2
    export kubo_oscillator_sde_3, kubo_oscillator_psde_3, kubo_oscillator_spsde_3
    export kubo_oscillator_ode

    q_init_A=[0.5, 0.0]
    q_init_B=[0.5 0.0 -0.5; 0.0 0.5 0.0]

    const noise_intensity = 0.1

    const Δt = 0.01
    const nt = 10


    function kubo_oscillator_sde_v(t, q, v_out)
        v_out[1]=  q[2]
        v_out[2]= -q[1]
    end


    function kubo_oscillator_sde_B(t, q, B_out::AbstractVector, ν=noise_intensity)
        B_out[1] = +ν*q[2]
        B_out[2] = -ν*q[1]
    end

    function kubo_oscillator_sde_B(t, q, B_out::AbstractMatrix, col=1, ν=noise_intensity)
        B_out[1,col] = +ν*q[2]
        B_out[2,col] = -ν*q[1]
    end


    function kubo_oscillator_sde_1(q₀=q_init_A)
        # q_init_A - single deterministic initial condition
        # Generating 3 sample paths
        # 1-dimensional noise
        SDE(1, 3, kubo_oscillator_sde_v, kubo_oscillator_sde_B, q₀)
    end

    function kubo_oscillator_sde_2(q₀=q_init_B)
        # q_init_B - interpreted as three random initial conditions
        # The 3 columns correspond to 3 sample paths
        # 1-dimensional noise
        SDE(1, 1, kubo_oscillator_sde_v, kubo_oscillator_sde_B, q₀)
    end

    function kubo_oscillator_sde_3(q₀=q_init_A)
        # q_init_A - interpreted as one random initial conditions with one sample path
        # 1-dimensional noise
        SDE(1, 1, kubo_oscillator_sde_v, kubo_oscillator_sde_B, q₀)
    end


    # ODE

    function kubo_oscillator_ode(q₀=q_init_A)
        ODE(kubo_oscillator_sde_v, q₀)
    end


    # PSDE

    q_init_C=[0.5]
    p_init_C=[0.0]

    q_init_D=[0.5 0.0 -0.5]
    p_init_D=[0.0 0.5 0.0]


    function kubo_oscillator_psde_v(t, q, p, v_out)
        v_out[1] =  p[1]
    end

    function kubo_oscillator_psde_f(t, q, p, f_out)
        f_out[1] = -q[1]
    end

    function kubo_oscillator_psde_B(t, q, p, B_out, ν=noise_intensity)
        B_out[1,1] = +ν*p[1]
    end

    function kubo_oscillator_psde_G(t, q, p, G_out, ν=noise_intensity)
        G_out[1,1] = -ν*q[1]
    end


    function kubo_oscillator_psde_1(q₀=q_init_C, p₀=p_init_C)
        # q_init_C - single deterministic initial condition
        # Generating 3 sample paths
        # 1-dimensional noise
        PSDE(1, 3, kubo_oscillator_psde_v, kubo_oscillator_psde_f, kubo_oscillator_psde_B, kubo_oscillator_psde_G, q₀, p₀)
    end

    function kubo_oscillator_psde_2(q₀=q_init_D, p₀=p_init_D)
        # q_init_D - interpreted as a single random initial condition
        # The 3 columns correspond to 3 sample paths
        # 1-dimensional noise
        PSDE(1, 1, kubo_oscillator_psde_v, kubo_oscillator_psde_f, kubo_oscillator_psde_B, kubo_oscillator_psde_G, q₀, p₀)
    end

    function kubo_oscillator_psde_3(q₀=q_init_C, p₀=p_init_C)
        # q_init_C - interpreted as a single random initial condition with one sample path
        # 1-dimensional noise
        PSDE(1, 1, kubo_oscillator_psde_v, kubo_oscillator_psde_f, kubo_oscillator_psde_B, kubo_oscillator_psde_G, q₀, p₀)
    end


    # SPSDE

    function kubo_oscillator_spsde_v(t, q, p, v_out)
        v_out[1] =  p[1]
    end

    function kubo_oscillator_spsde_f1(t, q, p, f_out)
        f_out[1] = -q[1]
    end

    function kubo_oscillator_spsde_f2(t, q, p, f_out)
        f_out[1] = 0
    end

    function kubo_oscillator_spsde_B(t, q, p, B_out, ν=noise_intensity)
        B_out[1,1] = +ν*p[1]
    end

    function kubo_oscillator_spsde_G1(t, q, p, G_out, ν=noise_intensity)
        G_out[1,1] = -ν*q[1]
    end

    function kubo_oscillator_spsde_G2(t, q, p, G_out, ν=noise_intensity)
        G_out[1,1] = 0
    end


    function kubo_oscillator_spsde_1(q₀=q_init_C, p₀=p_init_C)
        # q_init_C - single deterministic initial condition
        # Generating 3 sample paths
        # 1-dimensional noise
        SPSDE(1, 3, kubo_oscillator_spsde_v, kubo_oscillator_spsde_f1, kubo_oscillator_spsde_f2, kubo_oscillator_spsde_B, kubo_oscillator_spsde_G1, kubo_oscillator_spsde_G2, q₀, p₀)
    end

    function kubo_oscillator_spsde_2(q₀=q_init_D, p₀=p_init_D)
        # q_init_D - interpreted as a single random initial condition
        # The 3 columns correspond to 3 sample paths
        # 1-dimensional noise
        SPSDE(1, 1, kubo_oscillator_spsde_v, kubo_oscillator_spsde_f1, kubo_oscillator_spsde_f2, kubo_oscillator_spsde_B, kubo_oscillator_spsde_G1, kubo_oscillator_spsde_G2, q₀, p₀)
    end

    function kubo_oscillator_spsde_3(q₀=q_init_C, p₀=p_init_C)
        # q_init_C - interpreted as a single random initial condition with one sample path
        # 1-dimensional noise
        SPSDE(1, 1, kubo_oscillator_spsde_v, kubo_oscillator_spsde_f1, kubo_oscillator_spsde_f2, kubo_oscillator_spsde_B, kubo_oscillator_spsde_G1, kubo_oscillator_spsde_G2, q₀, p₀)
    end

end
