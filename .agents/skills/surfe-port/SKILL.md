---
name: surfe-port
description: Continue the fixed, dependency-ordered, one-task-per-session migration of frozen Surfe non-visual core into pure-Rust GeoRBF. Invoke explicitly with "$surfe-port continue"; never use for unrelated work.
---

# Surfe Port

将冻结在提交 `290dbe0ab344f4258a4935f05cad0f153f0f69a4` 的 Surfe 非可视化核心，按固定依赖顺序、一次一个任务地忠实迁移到 GeoRBF。仅在用户显式调用 `$surfe-port continue` 时执行本工作流。

## 权威来源

- 将仓库根目录 `docs/port/STATE.json` 作为跨会话任务状态的唯一事实来源。
- 将 `docs/port/PLAN.md` 作为固定任务定义、任务顺序和全局验收规则的权威来源。
- 将 `docs/port/JOURNAL.md` 作为已执行工作和证据的追加式日志。
- 不依赖旧聊天内容推断任务状态，不动态生成 Issue 或顶级任务。

## 1. 恢复上下文

每次调用 `$surfe-port continue` 时，按以下顺序执行：

1. 使用 `git rev-parse --show-toplevel` 定位仓库根目录，并从根目录工作。
2. 读取根目录 `AGENTS.md`；若文件不存在，记录这一事实后继续。
3. 完整读取：
   - `docs/port/STATE.json`
   - `docs/port/PLAN.md`
   - `docs/port/JOURNAL.md`
4. 只以 `STATE.json` 恢复任务状态，不依赖旧聊天。
5. 校验 `task_order`、`completed_tasks`、`last_completed_task`、`active_task` 和 `next_task` 一致：
   - `completed_tasks` 必须是 `INIT` 加 `task_order` 的连续前缀；
   - `last_completed_task` 必须等于 `completed_tasks` 的最后一项；
   - 正常待执行状态下 `active_task` 必须为 `null`，`next_task` 必须是该前缀后的唯一任务；
   - 中断恢复状态下 `active_task` 必须等于尚未完成的 `next_task`；
   - `status: complete` 仅允许在 T34 已完成且 `next_task` 为 `null` 时出现。
6. 从 `PLAN.md` 读取且只执行当前任务的正式定义，并读取完成该任务所必需的全局规则。不得提前执行后续任务。
7. 检查现有 GeoRBF 实现、测试和文档后再修改；不得假定仓库为空，不得无证据推倒重写。
8. 若 reference 是当前任务所需，按 `SURFE_REFERENCE_DIR`、`../surfe`、`.cache/surfe-reference` 的优先级解析，并校验其 HEAD 或检出状态为固定提交。无法获得 reference 时，不得伪造相关验证通过。

## 2. Git 安全

- 使用当前分支。
- 不创建分支、worktree、Issue 或 PR。
- 不 push。
- 不修改 Git 全局配置。
- 不使用 `git reset --hard`、`git clean` 或自动 stash。
- 开始任务前记录 `git status --short`，区分既有修改与本任务修改。
- 不修改或提交任务开始前存在的无关改动。
- 若存在无法与当前任务安全区分的 tracked 修改，停止并报告，不得擅自提交。
- 若工作区修改明显属于上次被中断的同一任务，则继续该任务，不得前进。
- 只 stage 当前任务相关文件以及同一任务更新的 `STATE.json`、`JOURNAL.md`。

## 3. 单任务执行

1. 若 `active_task` 为 `null`，将其设为 `next_task`；若已等于 `next_task`，恢复同一任务。任何时候不得选择其他任务。
2. 本会话只能处理该当前任务。复杂任务可以跨会话继续，但所有会话仍停留在同一任务，直至完成门槛全部通过。
3. 不提前实现下一任务，不进行顺手重构，不创建新的顶级任务，不更改固定顺序或扩大范围。
4. 对照固定 Surfe 提交阅读当前任务涉及的全部源码；在当前任务文档或日志中记录准确的 `path@290dbe0ab344f4258a4935f05cad0f153f0f69a4` 和相关符号。
5. 新发现若属于当前任务，在当前任务内解决；若属于后续内容，只登记到对应的现有后续任务，不创建新任务。
6. 先增加或更新测试，再完成实现。
7. 遵守纯 Rust 边界：生产构建、测试和发布不得依赖 C++、Eigen、Qt、VTK、CMake、bindgen、`cxx`、生产用 `cc`、MKL、OpenBLAS、LAPACK 或 native BLAS FFI。
8. 排除 Qt、VTK、GUI、`geo_builder`、等值面显示和纯可视化代码。冻结 C++ 只可存在于仓库外部或被忽略的 reference/oracle 目录，不得进入生产源码或发布包。
9. 忠实保留有效输入的约束语义、核与导数、矩阵布局、求解路径、标量场、梯度场、层位关系、错误分类和可达功能。不得复制 C++ 的未定义行为、内存问题或数据竞争。
10. 发现源码缺陷时，先用兼容测试记录可观察行为，再在既有兼容文档中明确处理；不得静默改变数学结果。
11. 不得仅因条件数很差而提前拒绝系统；按冻结 Surfe 的实际语义尝试求解，再依据有限性、残差和约束可行性判断。
12. 全局 parity 完成前，不做改变数学语义的算法重构，不引入 FMM、低秩、局部近似、稀疏替代或新的求解定义。
13. 运行当前任务在 `PLAN.md` 定义的全部验证。测试未执行、reference 无法运行、parity 未证实或性能未判断时，不得写成通过。
14. 任一完成门槛失败时，保持 `active_task` 和 `next_task` 为当前任务，不得更新完成列表，不得提交“完成”提交。

## 4. 完成与提交

只有当前任务全部验收通过后，才执行以下操作：

1. 更新 `STATE.json`：
   - 将当前任务追加到 `completed_tasks`；
   - 将 `last_completed_task` 更新为当前任务；
   - 将 `active_task` 设为 `null`；
   - 将 `next_task` 设为固定序列中的唯一后继；
   - T34 完成后将 `next_task` 设为 `null`，将 `status` 设为 `complete`，并仅依据实测证据更新 `final_definition_of_done`。
2. 向 `JOURNAL.md` 追加：
   - 任务编号和名称；
   - 阅读的 Surfe 源码及准确提交；
   - 修改文件；
   - 核心实现；
   - 验证命令和结果；
   - parity 证据；
   - 性能证据；
   - 已记录但属于现有后续任务的发现。
3. 只 stage 当前任务相关文件、`STATE.json` 和 `JOURNAL.md`。
4. 自动提交，使用格式 `port(surfe): complete TNN <concise-task-name>`。INIT 专用格式为 `chore(port): initialize surfe-port controller`。
5. 检查提交成功，取得 short SHA 和 subject，并确认没有遗留由本任务产生但未提交的文件。
6. 不 push。

## 5. 输出格式

普通任务完成后只输出以下五段，最后一行必须恰好是 `$surfe-port continue`：

```text
完成：TNN — <任务名称>
提交：<short-sha> <commit subject>
验证：<最重要的验证结果>
下一任务：TNN — <唯一后继任务名称>

$surfe-port continue
```

只能给出一个下一任务，不提供多个选项，不要求用户复制其他长提示词，不主动建议并行任务。

任务未完成时只输出：

```text
未完成：TNN — <任务名称>
阻塞：<准确原因>
状态：仍停留在 TNN

$surfe-port continue
```

不得推进状态。

T34 完成后只输出：

```text
完成：T34 — 最终发布验收
提交：<short-sha> <commit subject>
状态：COMPLETE
结论：GeoRBF 已达到冻结 Surfe 非可视化核心行为等价、纯 Rust 和性能门槛。
```

最终状态下不得再输出下一任务命令。
