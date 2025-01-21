\documentclass{article}
\usepackage{amsmath}
\usepackage{amssymb}
\usepackage{graphicx}
\usepackage{enumitem}

\title{Performance Metrics for Sparse Blind Deconvolution}
\author{Dong Lab}
\date{\today}

\begin{document}
\maketitle

\section{Performance Metrics}
To evaluate the performance of our algorithm, we employ four distinct metrics:

\begin{enumerate}[label=\arabic*.]
    \item \textbf{Kernel Quality:} This metric quantifies the accuracy of kernel recovery by computing the normalized inner product between the recovered and ground truth kernels:
    \begin{equation}
        Q_k = \frac{|\langle A, A_0 \rangle|^2}{\|A\|^2_F \|A_0\|^2_F}
    \end{equation}
    where $A$ is the recovered kernel, $A_0$ is the ground truth kernel, and $\|\cdot\|_F$ denotes the Frobenius norm. This metric ranges from 0 to 1, with 1 indicating perfect recovery.

    \item \textbf{Activation Recovery:} This measures the accuracy of activation map reconstruction after accounting for potential spatial shifts:
    \begin{equation}
        Q_a = \frac{|\langle X_{\text{aligned}}, X_0 \rangle|^2}{\|X_{\text{aligned}}\|^2_F \|X_0\|^2_F}
    \end{equation}
    where $X_{\text{aligned}}$ represents the recovered activation map after optimal alignment with the ground truth $X_0$. Values range from 0 to 1, with 1 indicating perfect recovery.

    \item \textbf{Demixing Score:} This quantifies how well the algorithm separates different components by measuring the normalized cross-correlation between different activation maps:
    \begin{equation}
        D = 1 - \frac{1}{N(N-1)} \sum_{i\neq j} \max_{\tau} \left|\frac{\langle X_i, X_j(\tau) \rangle}{\|X_i\|_F \|X_j\|_F}\right|
    \end{equation}
    where $N$ is the number of kernels, $X_i$ and $X_j$ are recovered activation maps, and $\tau$ represents all possible spatial shifts. A score of 1 indicates perfect separation (no cross-correlation), while 0 indicates complete mixing.

    \item \textbf{Runtime:} The computational efficiency is measured by the total execution time:
    \begin{equation}
        T = t_{\text{end}} - t_{\text{start}}
    \end{equation}
    measured in seconds, including all iterations and phases of the algorithm.
\end{enumerate}

\section{Interpretation}
These metrics provide complementary measures of algorithm performance:
\begin{itemize}
    \item The first three metrics ($Q_k$, $Q_a$, and $D$) are normalized to [0,1], where higher values indicate better performance
    \item Runtime ($T$) is measured in absolute terms (seconds), with lower values being preferable
    \item Kernel Quality and Activation Recovery require ground truth for computation
    \item The Demixing Score can be computed without ground truth, making it suitable for real-world applications
\end{itemize}

\section{Implementation Notes}
The metrics are implemented in MATLAB with the following considerations:
\begin{itemize}
    \item Cross-correlations are computed using \texttt{normxcorr2} for efficiency
    \item Alignment for activation recovery uses maximum correlation to determine optimal shifts
    \item Runtime includes all preprocessing, iterations, and postprocessing steps
    \item All metrics are computed per kernel and then averaged for multi-kernel scenarios
\end{itemize}

\end{document}