# 项目协作规范

## 🌿 分支管理
- `main`：稳定版本（仅管理员可合并）
- `dev`：开发分支（主要工作分支）
- `feat/*`：功能开发分支（如 `feat/variant-annotation`）

## 💡 代码提交
1. 新建功能分支：  
   `git checkout -b feat/your-feature`
2. 提交清晰的信息：  
   `git commit -m "fix: correct hisat2 parameters in pipeline"`
3. 推送到远程仓库：  
   `git push origin feat/your-feature`
4. 创建Pull Request到`dev`分支

## ✅ 代码审查
- 至少需要1名协作者审查
- 通过CI测试后才能合并
